# [[file:../fresa.org::testrosenbrock][testrosenbrock]]
using Fresa
nx = 2
x0 = 0.5*ones(nx)
# specify the domain for the search, x ∈ [0,10]ⁿ
d = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))
rosenbrock(x) = ([(1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2], 0)
# f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x0,rosenbrock)]
# now invoke Fresa to solve the problem
best, pop = Fresa.solve(rosenbrock, p0;
                        domain=d,
                        ngen=1000,
                        issimilar = Fresa.similarx,
                        ϵ=1e-8)
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
println("identified after $(best.since[1]) function evaluates in generation $(best.since[2]).")
# testrosenbrock ends here
