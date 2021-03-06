# [[file:../fresa.org::testsimple][testsimple]]
# load in the Fresa optimization package
using Fresa
# specify the dimension of the search space
nx = 2
# create an initial point in the search space
x0 = 0.5*ones(nx)
# specify the domain for the search, x ∈ [0,10]ⁿ
domain = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))
# the actual objective function
f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0.0)
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x0,f)]
# now invoke Fresa to solve the problem
@time best, pop = Fresa.solve(f, p0, domain)
# output the results
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testsimple ends here

# [[file:../fresa.org::*simple objective function][simple objective function:2]]
#+plot: ind:1 deps:(2) with:"linespoints lt 3 pt 7 ps 0.25" set:nokey set:"yrange [0:1]" set:"xrange [0:*]"
println("#+plot: ind:1 deps:(2) with:\"linespoints pt 7 ps 0.25\" set:nokey set:\"yrange [0:1]\" set:\"xrange [0:*]\" set:\"xlabel 'Generation'\" set:\"ylabel 'fitness'\"")
Fresa.printHistoryTrace(best)
# simple objective function:2 ends here
