using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JuMP, Ipopt, Juniper, Cbc

# <editor-fold> Initial things
### Parameters
global F = 1:1 # Set of followers
global L = 1:3 # set of leaders
lenF = length(F)
lenL = length(L)
α = 10
c_cst = 1 # Marginal costs are the same for everyone
MC_F = fill(c_cst, lenF)
MC_L = fill(c_cst, lenL)
M = 1e3 # Big M parameter
@enum DiagonalisationType GaussSeidel Jacobi
diagType = GaussSeidel
# diagType = Jacobi
# NOTE: Jacobi method seems to be unstable???
# </editor-fold>

# <editor-fold> Create initial model
### Define Model and associated Optimizer
# MINLP, so slightly weird parameters for optimizer
optimizer = Juniper.Optimizer
params = Dict{Symbol,Any}()
params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level = 0)
params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel = 0)
m = Model(with_optimizer(optimizer, params))

### Variables
@variable(m, Q[L] >= 0, start = 5)
@variable(m, q[F] >= 0)
@variable(m, b[F], Bin)

### Constraints
# These are seen by all the leaders
@constraint(m, [f=F],
    MC_F[f] - α + 2*q[f]
    + sum(
        q[ff] for ff in F if ff != f
    )
    + sum(
        Q[l] for l in L
    )
    >= 0
)
@constraint(m, [f=F],
    q[f] <= M*(1 - b[f])
)
@constraint(m, [f=F],
    MC_F[f] - α + 2*q[f]
    + sum(
        q[ff] for ff in F if ff != f
    )
    + sum(
        Q[l] for l in L
    )
    <= M*b[f]
)
# </editor-fold>

# <editor-fold> Loop to solve
# Initialise Q values - rows = Q, columns = iteration
global storedQ = rand(lenL, 1)

# Loop until tolerance met or maximum iterations
tolerance = 0.01
maxIter = 10
while size(storedQ, 2) <= 1 ||
        (
            all(
                abs.(storedQ[:, end - 1] .- storedQ[:, end]) .> tolerance
            ) && size(storedQ, 2) <= maxIter
        )
    # Julia doesn't declare global variables in while loops (???)
    global L, F, storedQ

    # Repeat last row of storedQ
    storedQ = hcat(storedQ, storedQ[:, end])

    # Loop over the leaders and solve the optimisation problem each time
    for l in L
        # Unfix the leaders, but also make sure they're positive!
        for l in L
            if is_fixed(Q[l])
                unfix(Q[l])
                set_lower_bound(Q[l], 0)
            end
        end

        # Set objective
        @NLobjective(
            m, Max, (α -
            (
                sum(Q[l] for l in L) + sum(q[f] for f in F)
            )
            )*Q[l]- MC_L[l]*Q[l]
        )

        # Fix the values of all the leaders who aren't being considered
        MM = [ll for ll in L if ll != l]
        if diagType == GaussSeidel
            for ll in MM
                fix(Q[ll], storedQ[ll,end], force = true)
            end
        elseif diagType == Jacobi
            for ll in MM
                fix(Q[ll], storedQ[ll, end - 1], force = true)
            end
        else
            throw()
        end


        # Solve the optimisation problem
        optimize!(m)

        # Store the solved value of Q[l]
        storedQ[l,end] = value(Q[l])
    end
end
# </editor-fold>

# Show storedQ
println("""
Leaders = $(storedQ[:, end]))
Followers = $(value.(q).data)
""")

σ
