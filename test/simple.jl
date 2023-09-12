# [[file:../fresa.org::testsimple][testsimple]]
# load in the Fresa optimization package
using Fresa

# specify the dimension of the search space
nx = 2

# create an initial point in the search space
x0 = 0.5*ones(nx)

# specify the domain for the search, x ∈ [0,10]ⁿ, by specifying fixed
# lower and upper bounds.  This will test the creation of the Domain
# data type by the solve method
a = zeros(length(x0))
b = 10.0 * ones(length(x0))

# The alternatie would be to define the Domain directly:
# d = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))

# the actual objective function
f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0)

# create the initial population consisting of this single point
p0 = [Fresa.Point(x0,f)]

# now invoke Fresa to solve the problem
best, pop = Fresa.solve(f, p0, lower = a, upper = b)

# output the results
println("Population at end:\n$pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testsimple ends here

# [[file:../fresa.org::*simple objective function][simple objective function:2]]
println("\nHistory trace, by generation number, of fitness value of solution selected for propagation which results in a new best solution:")
println("#+plot: ind:1 deps:(2) with:\"linespoints pt 7 ps 0.25\" set:nokey set:\"yrange [0:1]\" set:\"xrange [0:*]\" set:\"xlabel 'Generation'\" set:\"ylabel 'fitness'\"")
Fresa.printHistoryTrace(best)
# simple objective function:2 ends here
