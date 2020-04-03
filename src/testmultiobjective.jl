# [[file:~/s/research/julia/Fresa.jl/src/fresa.org::testmultiobjective][testmultiobjective]]
using Fresa
nx = 2
a = zeros(nx)
b = ones(nx)
x = rand(nx)
f = x -> ( [sin(x[1]-x[2]); cos(x[1]+x[2])], 0)
pareto, population = Fresa.solve(f, x, a, b;
                                 #fitnesstype = :hadamard,
                                 #fitnesstype = :borda,
                                 fitnesstype = :nondominated,
                                 ngen=100,
                                 npop=10,
                                 plotvectors=true,
                                 tolerance=0.01)

println("Pareto front:")
println(population[pareto])
#using BenchmarkTools
#@benchmark

using PyPlot
z = [population[pareto[i]].z for i in 1:length(pareto)];
PyPlot.plot([z[i][1] for i=1:length(z)],
            [z[i][2] for i=1:length(z)],
            "ro")
PyPlot.savefig("x.pdf")
# testmultiobjective ends here
