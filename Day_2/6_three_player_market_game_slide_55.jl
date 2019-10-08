using Complementarity, JuMP

α = 10 # intercept of demand curve, highest price you would pay
β = 1 # slope of the demand curve
γ1 = 1 # Marginal cost of producers
γ2 = 1
γ3 = 1
N = 3

m = MCPModel()

@variable(m, q1 >= 0)
@variable(m, q2 >= 0)
@variable(m, q3 >= 0)

@mapping(
    m, F1, N*q1 - β*(q2 + q3) + γ1 - α
)
@mapping(
    m, F2, N*q2 - β*(q1 + q3) + γ2 - α
)
@mapping(
    m, F3, N*q3 - β*(q1 + q2) + γ3 - α
)

@complementarity(m, F1, q1)
@complementarity(m, F2, q2)
@complementarity(m, F3, q3)

solveLCP(m)

println("""
q1 = $(result_value(q1))
q2 = $(result_value(q2))
q3 = $(result_value(q3))
""")
