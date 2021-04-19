# [[file:../fresa.org::testmultiobjective][testmultiobjective]]
using Fresa
nx = 2
a = zeros(nx)
b = ones(nx)
x = rand(nx)
f = x -> ( [sin(x[1]-x[2]); cos(x[1]+x[2])], 0)
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x,f)]
# now invoke Fresa to solve the problem
pareto, population = Fresa.solve(f, p0, a, b;
                                 #fitnesstype = :hadamard,
                                 #fitnesstype = :borda,
                                 fitnesstype = :nondominated,
                                 ngen=200,
                                 npop=20,
                                 plotvectors=true,
                                 tolerance=0.01)

println("**** Pareto front:")
println("#+plot: ind:2 deps:(3) with:points")
println(population[pareto])
#using BenchmarkTools
#@benchmark
# testmultiobjective ends here
