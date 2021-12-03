# [[file:../fresa.org::testrosenbrock][testrosenbrock]]
using Fresa
nx = 2
x0 = 0.5*ones(nx)
# specify the domain for the search, x ∈ [0,10]ⁿ
domain = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))
rosenbrock(x) = ([(1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2], 0)
# f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0)
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x0,rosenbrock)]
# now invoke Fresa to solve the problem
best, pop = Fresa.solve(rosenbrock, p0, domain; ngen=1000, tolerance=1e-8)
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testrosenbrock ends here
