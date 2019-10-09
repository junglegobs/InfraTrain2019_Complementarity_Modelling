using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JuMP, Ipopt, Juniper, Cbc

### Parameters
global N = 1:1 # Set of followers
global M = 1:3 # set of leaders
lN = length(N)
lM = length(M)
α = 10
c_cst = 1
c = fill(c_cst, lN)
C = fill(c_cst, lM)
K = 1e3
@enum DiagonalisationType GaussSeidel Jacobi
diagType = GaussSeidel
# diagType = Jacobi
# NOTE: Jacobi method seems to be unstable???

### MINLP, so slightly weird parameters for optimizer
optimizer = Juniper.Optimizer
params = Dict{Symbol,Any}()
params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level = 0)
params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel = 0)
m = Model(with_optimizer(optimizer, params))

### Variables
@variable(m, Q[M] >= 0, start = 5)
@variable(m, q[N] >= 0)
@variable(m, b[N], Bin)

### Constraints
# These are seen by all the leaders
@constraint(m, [f=N],
    c[f] - α + 2*q[f]
    + sum(
        q[ff] for ff in N if ff != f
    )
    + sum(
        Q[l] for l in M
    )
    >= 0
)
@constraint(m, [f=N],
    q[f] <= K*(1 - b[f])
)
@constraint(m, [f=N],
    c[f] - α + 2*q[f]
    + sum(
        q[ff] for ff in N if ff != f
    )
    + sum(
        Q[l] for l in M
    )
    <= K*b[f]
)

# Initialise Q values - rows = Q, columns = iteration
global storedQ = rand(lM, 1)

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
    global M, N, storedQ

    # Repeat last row of storedQ
    storedQ = hcat(storedQ, storedQ[:, end])

    # Loop over the leaders and solve the optimisation problem each time
    for l in M
        # Unfix the leaders, but also make sure they're positive!
        for l in M
            if is_fixed(Q[l])
                unfix(Q[l])
                set_lower_bound(Q[l], 0)
            end
        end

        # Set objective
        @NLobjective(
            m, Max, (α -
            (
                sum(Q[l] for l in M) + sum(q[f] for f in N)
            )
            )*Q[l]- C[l]*Q[l]
        )

        # Fix the values of all the leaders who aren't being considered
        MM = [ll for ll in M if ll != l]
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

# Show storedQ
println("""
Leaders = $(storedQ[:, end]))
Followers = $(value.(q).data)
""")
