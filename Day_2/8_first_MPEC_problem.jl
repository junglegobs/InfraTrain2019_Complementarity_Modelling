using JuMP, Ipopt, Juniper, Cbc

N = 1
α = 10
c_followers = 1
c = fill(c_followers, N)
C = 1
K = 1e3

# MINLP, so slightly weird parameters
# https://github.com/lanl-ansi/Juniper.jl
optimizer = Juniper.Optimizer
params = Dict{Symbol,Any}()
params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level = 0)
params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel = 0)
m = Model(with_optimizer(optimizer, params))

# Variables
@variable(m, Q >= 0, start = 5)
@variable(m, q[1:N] >= 0)
@variable(m, b[1:N], Bin)

# Constraints
@constraint(m, [i=1:N],
    c[i] - α + 2*q[i] + sum(
        q[j] for j = 1:N if j != i
    ) + Q >= 0
)
@constraint(m, [i=1:N],
    q[i] <= K*(1 - b[i])
)
@constraint(m, [i=1:N],
    c[i] - α + 2*q[i] + sum(
        q[j] for j = 1:N if j != i
    ) + Q <= K*b[i]
)

# Objective
@NLobjective(m, Max,
    (α - (
            Q + sum(q[i] for i=1:N)
        )
    )*Q - C*Q
)

# Solve!
optimize!(m)

# Print the solution
q_total = value(Q) + sum(value.(q))
price = α - sum(q_total)
println("""Leader production = $(value(Q))""")
println("""Total production = $q_total""")
println("""Price = $price""")
