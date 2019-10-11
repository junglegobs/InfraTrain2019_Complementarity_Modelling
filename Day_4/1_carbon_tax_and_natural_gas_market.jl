using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JuMP, Ipopt, Plots

# <editor-fold> Initial things
### Parameters
len_P = 10
len_NG = 7
len_BG = 3
P = 1:len_P # Set of producers
NG = 1:len_NG # Set of natural gas
BG = len_NG + 1:len_BG # set of biogas producers

# Emission factors
EF_BG = 0.0
EF_NG = 0.48
EF = [p in NG ? EF_NG : EF_BG for p in P]
EF_bin = [p in NG ? true : false for p in P]

# Production marginal costs
c_BG = 0.7
c_NG = 0.5
c = [p in NG ? c_NG : c_BG for p in P]

# Capacities
cap = [100 for p in P]

# Other parameters
α = 10

# </editor-fold>

# <editor-fold> Define model
m = Model(with_optimizer(Ipopt.Optimizer))

# Define variables
@variable(m, q[P] >= 0) # production
@variable(m, t >= 0) # tax
@variable(m, λ[P] >= 0) # Dual of capacity constraint

# Expressions
total_emissions = @NLexpression(m,
    sum(
        EF[p]*q[p] for p in P
    )
)
producer_profit = @NLexpression(m,
    sum(
        - c[p]*q[p] + (α - sum(q[p] for p in P))*q[p] - t*q[p]*EF_bin[p] for p = P
    )
)
consumer_surplus = @NLexpression(m,
    0.5 * sum(q[p] for p in P)^2
)
tax_revenue = @NLexpression(m,
    sum(t*q[p]*EF_bin[p] for p = P)
)
social_welfare = @NLexpression(m,
    producer_profit + consumer_surplus + tax_revenue
)

# Inequality constraints
@constraint(m, [p=P],
    c[p] - α + 2*q[p] + sum(q[pp] for pp in P if pp != p)+ λ[p] + t*EF_bin[p] >= 0
)
@constraint(m, [p=P],
    q[p] <= cap[p]
)

# Inner product constraints
@NLconstraint(m, [p=P],
    (
        c[p] - α + 2*q[p] + sum(q[pp] for pp in P if pp != p)+ λ[p] + t*EF_bin[p]
    ) * (
        q[p]
    )
    == 0
)
@NLconstraint(m, [p=P],
    (cap[p] - q[p])*(λ[p]) == 0
)
# </editor-fold>

# <editor-fold> Weighted method loop for multiobjective
# w = 0, min total emissions, w = 1, max social welfare
obj_weights = [0, 1]
tax_store = fill(0.0, length(obj_weights))
total_emissions_store = fill(0.0, length(obj_weights))
social_welfare_store = fill(0.0, length(obj_weights))
total_production_store = fill(0.0, length(obj_weights))

for i in 1:length(obj_weights)
    w = obj_weights[i]

    # formulate objective
    println("\n\n\n $w \n\n\n")
    @NLobjective(m, Max, w*social_welfare - (1 - w)*total_emissions)

    # solve model
    optimize!(m)

    # store tax values
    tax_store[i] = value(t)
    total_production_store[i] = sum(value.(q).data)
    total_emissions_store[i] = value(total_emissions)
    social_welfare_store[i] = value(social_welfare)
end
# </editor-fold>

# <editor-fold> Weighted method loop for multiobjective
max_allowed_emissions_range = range(
    minimum(total_emissions_store),
    maximum(total_emissions_store),
    length = 40
)
tax_store = fill(0.0, length(max_allowed_emissions_range))
total_emissions_store = fill(0.0, length(max_allowed_emissions_range))
social_welfare_store = fill(0.0, length(max_allowed_emissions_range))
total_production_store = fill(0.0, length(max_allowed_emissions_range))

# Add emissions constraint
@NLparameter(m, max_allowed_emissions == max_allowed_emissions_range[1])
@NLconstraint(m,
    total_emissions <= max_allowed_emissions
)
# formulate objective
@NLobjective(m, Max, social_welfare)

for i in 1:length(max_allowed_emissions_range)
    # Fix maximum emissions
    set_value(max_allowed_emissions, max_allowed_emissions_range[i])

    # solve model
    optimize!(m)

    # store tax values
    tax_store[i] = value(t)
    total_production_store[i] = sum(value.(q).data)
    total_emissions_store[i] = value(total_emissions)
    social_welfare_store[i] = value(social_welfare)
end
# </editor-fold>

# <editor-fold> Plot tax

plot(social_welfare_store, total_emissions_store,
    xlabel = "Social welfare",
    ylabel = "Total emissions",
    seriestype = :scatter,
    lab = "",
    xlims = (39,44.5)
)

# </editor-fold>
