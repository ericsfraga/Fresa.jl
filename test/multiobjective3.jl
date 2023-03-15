# [[file:../fresa.org::testmultiobjective3][testmultiobjective3]]
using Fresa
using Profile
nx = 5
# specify the domain for the search, x ∈ [0,1]ⁿ
d = Fresa.Domain(x -> zeros(length(x)), x -> ones(length(x)))
x = zeros(nx)
f = x -> ([ sum((x.-0.5).^2 .+ 1)
            sum(cos.(x))
            sum(sin.(x))],
          0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x,f)]
# now invoke Fresa to solve the problem
pareto, population = Fresa.solve(f, p0;
                                 domain = d,
                                 ϵ = 0.01,
                                 archiveelite = false,
                                 issimilar = Fresa.similarx,
                                 np=20,
                                 ngen=300
                                 )

println("*** Pareto front")
println("Pareto set of size $(length(pareto)) with indices: $pareto")
println(population[pareto])
# testmultiobjective3 ends here
