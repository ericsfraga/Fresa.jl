# [[file:../fresa.org::testsimple][testsimple]]
using Fresa
nx = 2
x0 = 0.5*ones(nx)
a = zeros(nx)
b = 10*ones(nx)
f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0.0)
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x0,f)]
# now invoke Fresa to solve the problem
@time best, pop = Fresa.solve(f, p0, a, b)
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testsimple ends here
