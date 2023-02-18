# [[file:fresa.org::examplesolution][examplesolution]]
using Fresa
    function objective(x)
        # calculate the objective function value
        z = 5*x[1]^2 + 4*x[2]^2 - 60*x[1] - 80*x[2]
        # evaluate the constraints so that feasible points result in a
        # non-positive value, i.e. 0 or less, but infeasible points give a
        # positive value.  We choose the maximum of both constraints as
        # the value to return as an indication of feasibility
        g = maximum( [ 6*x[1] + 5*x[2] - 60
                       10*x[1] + 12*x[2] - 150 ] )
        # return the objective function value along with indication of
        # feasibility
        (z, g)
    end
    dom = Fresa.Domain( x -> [ 0.0,  0.0 ],  # lower bounds
                        x -> [ 8.0, 12.5 ] ) # upper bounds
    initialpopulation = [ Fresa.Point( [4.0, 6.25 ], objective ) ]
    best, population = Fresa.solve( objective, # function 
                                    initialpopulation, # initial points
                                    domain = dom )     # the search domain
    println("Population at end:")
    println("$population")
    println("Best solution found is:")
    println("  f($( best.x ))=$( best.z )")
    println("with constraint satisfaction (≤ 0) or violation (> 0):")
    println("  g=$( best.g ).")
# examplesolution ends here
