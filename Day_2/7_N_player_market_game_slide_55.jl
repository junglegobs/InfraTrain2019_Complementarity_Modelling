using Complementarity, JuMP

N = 60 # Number of players
α = 10 # intercept of demand curve, highest price you would pay
β = 1 # slope of the demand curve
c = 1 # Marginal cost
γ = fill(c, N) # Marginal cost of producers

# Make model
m = MCPModel()

# Make production variables
@variable(m, q[i=1:N] >= 0)

# Make mapping of function F
# Remember: 0 ≤ F(x) ⟂ x ≥ 0
# In this case, x => q
# And F[i] => γ[i] - α + 2*q[i] + β*(∑ q[j] for j = 1:N, j != i)
@mapping(
    m, F[i=1:N],
    γ[i] - α + 2*q[i] + β*(
        sum(
            q[j] for j = 1:N if j != i
        )
    )
)

# Define complementarity
@complementarity(m, F, q)

# solve (need to use NLSolve for problems greater than 40 people)
solveLCP(m, solver=:NLsolve)

# Calculate some additional results
q_total = sum(
    result_value.(q)
)
price = α - β*sum(
    result_value.(q)
)
profit = price*result_value(q[1]) - c*result_value(q[1])

# display results
for i in 1:N
    println("""
    q$i = $(result_value(q[i]))
    """)
end
println("""Total production: $q_total""")
println("""Price: $price""")
println("""Profit: $profit""")
