using JuMP, Gurobi, Ipopt

m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m, x, start = pi/2) # hot start
@constraint(m, -x <= 0)
@constraint(m, x <= 4*pi)
@NLobjective(m, Min, sin(x))

optimize!(m)
println("""\n
x = $(value(x)/pi)/π
objective = $(objective_value(m))
""")

set_start_value(x, pi)
