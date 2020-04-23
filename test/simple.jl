# [[file:~/s/research/julia/Fresa.jl/fresa.org::testsimple][testsimple]]
using Fresa
nx = 2
x0 = 0.5*ones(nx)
a = zeros(nx)
b = 10*ones(nx)
f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0.0)
@time best, pop = Fresa.solve(f, x0, a, b)
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
# testsimple ends here
