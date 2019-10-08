using JuMP, Ipopt

m = Model(with_optimizer(Ipopt.Optimizer))

# Parameters
α = 10
c1 = 1
c2 = 1

# Variables
@variable(m, q1)
@variable(m, q2)

# Objective
@objective(m, Max,
    (α - (q1 + q2))*q1 - c*q1
    + (α - (q1 + q2))*q2 - c2*q2
)

q_start_values = [
    [1 1],
    [5 5],
    [0 0]
]

for q in q_start_values
    set_start_value(q1, q[1])
    set_start_value(q2, q[2])

    # Optimize
    optimize!(m)

    println("""
    q1 = $(value(q1))
    q2 = $(value(q2))
    p = $(α - value(q1))
    p = $(α - value(q2))
    """)
end

# Always get the same value, but actually q1 and q2 can be anything as long as q1 + q2
