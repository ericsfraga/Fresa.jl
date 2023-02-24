# [[file:../fresa.org::testmultiobjective][testmultiobjective]]
using Fresa
nx = 2
# specify the domain for the search, x ∈ [0,10]ⁿ
d = Fresa.Domain(x -> zeros(length(x)), x -> ones(length(x)))
# initial point in domain
x = rand(nx)
# objective function 
f = x -> ( [sin(x[1]-x[2]); cos(x[1]+x[2])], 0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x,f)]
# now invoke Fresa to solve the problem
pareto, population = Fresa.solve(f, p0;
                                 domain = d,
                                 ϵ = 0.01,
                                 #fitnesstype = :hadamard,
                                 #fitnesstype = :borda,
                                 fitnesstype = :nondominated,
                                 issimilar = Fresa.similarx,
                                 ngen=200,
                                 np=(20,40),
                                 plotvectors=true)

println("**** Pareto front:")
println("#+plot: ind:1 deps:(2) with:points")
println(population[pareto])
#using BenchmarkTools
#@benchmark
# testmultiobjective ends here
