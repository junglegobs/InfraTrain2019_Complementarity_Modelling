using JuMP, Ipopt

m = Model(with_optimizer(Ipopt.Optimizer))

# Parameters
α = 10
c = 1

@variable(m, q)

@objective(m, Max, (α - q)*q - c*q)

optimize!(m)

println("""
q = $(value(q))
""")
