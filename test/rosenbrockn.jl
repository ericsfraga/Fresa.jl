# [[file:../fresa.org::testrosenbrockn][testrosenbrockn]]
using Fresa
nx = 20
x0 = 0.5*ones(nx)
# specify the domain for the search, x ∈ [0,10]ⁿ
d = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))
rosenbrock(x) = (sum([100 * (x[i+1]-x[i]^2)^2 + (1-x[i])^2 for i ∈ 1:length(x)-1]), 0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x0,rosenbrock)]
# now invoke Fresa to solve the problem
best, pop = Fresa.solve(rosenbrock, p0;
                        domain=d,
                        np=100,
                        ngen=1000,
                        ϵ=1e-8,
                        issimilar = Fresa.similarx,
                        multithreading=true,
                        ticker = false) # test this option
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testrosenbrockn ends here
