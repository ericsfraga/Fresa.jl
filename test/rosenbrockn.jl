# [[file:../fresa.org::testrosenbrockn][testrosenbrockn]]
using Fresa
nx = 20
x0 = 0.5*ones(nx)
# specify the domain for the search, x ∈ [0,10]ⁿ
domain = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))
rosenbrock(x) = (sum([100 * (x[i+1]-x[i]^2)^2 + (1-x[i])^2 for i ∈ 1:length(x)-1]), 0)
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x0,rosenbrock)]
# now invoke Fresa to solve the problem
best, pop = Fresa.solve(rosenbrock, p0, domain; npop=100, ngen=1000, tolerance=1e-8, multithreading=true)
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testrosenbrockn ends here