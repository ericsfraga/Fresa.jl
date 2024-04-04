# [[file:fresa.org::examplesolution][examplesolution]]
using Fresa
function objective(x)
    # calculate the objective function value
    z = 5*x[1]^2 + 4*x[2]^2 - 60*x[1] - 80*x[2]
    # evaluate the constraints so that feasible points result in a
    # non-positive value, i.e. 0 or less, but infeasible points give a
    # positive value.  We choose the maximum of both constraints as
    # the value to return as an indication of feasibility
    g = max( 6*x[1] + 5*x[2] - 60,
             10*x[1] + 12*x[2] - 150 )
    # return the objective function value along with indication of
    # feasibility
    (z, g)
end
a = [ 0.0, 0.0 ]
b = [ 8.0, 12.5 ]
initialpopulation = [ Fresa.Point( [4.0, 6.25 ], objective ) ]
best, population = Fresa.solve( objective,         # the function 
                                initialpopulation, # initial points
                                lower = a,         # lower bounds
                                upper = b          # upper bounds
                                )
println("Population at end:")
println("$population")
println("Best solution found is:")
println("  f($( best.x ))=$( best.z )")
println("with constraint satisfaction (≤ 0) or violation (> 0):")
println("  g=$( best.g ).")
# examplesolution ends here
