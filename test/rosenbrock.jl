# [[file:~/s/research/julia/Fresa.jl/fresa.org::testrosenbrock][testrosenbrock]]
using Fresa
nx = 2
x0 = 0.5*ones(nx)
a = zeros(nx)
b = 10*ones(nx)
rosenbrock(x) = ([(1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2], 0)
# f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0)
best, pop = Fresa.solve(rosenbrock, x0, a, b)
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testrosenbrock ends here
