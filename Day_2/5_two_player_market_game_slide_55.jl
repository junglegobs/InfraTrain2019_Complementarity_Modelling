using Complementarity, JuMP

α = 10
Β = 1
γ1 = 1
γ2 = 1

m = MCPModel()

@variable(m, q1 >= 0)
@variable(m, q2 >= 0)

@mapping(m, F1, γ1 - (α - 2*q1 -q2))
@mapping(m, F2, γ2 - (α - 2*q2 - q1))

@complementarity(m, F1, q1)
@complementarity(m, F2, q2)

solveLCP(m)

println("""
q1 = $(result_value(q1))
q2 = $(result_value(q2))
""")
