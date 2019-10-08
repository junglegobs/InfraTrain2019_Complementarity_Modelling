# example from slide 46 version 2 (defining equations line by line)

# including packages
using Complementarity

# Same model in a more readable representation
m = MCPModel()

#adding variables
x1 = @variable(m, x1 >= 0)
x2 = @variable(m, x2 >= 0)
y1 = @variable(m, y1)

#setting up camplementary equations
F1 = @mapping(m, map1, x1 + x2)
F2 = @mapping(m, map2, x1 - y1)
F3 = @mapping(m, map3, x1 + x2 + y1 - 2)
@complementarity(m, map1, x1)
@complementarity(m, map2, x2)
@complementarity(m, map3, y1)

#soving model
solveMCP(m, linear=true)

@show result_value.(x1)
@show result_value.(x2)
@show result_value.(y1)
