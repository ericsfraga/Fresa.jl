# [[file:../fresa.org::testmultiobjective3][testmultiobjective3]]
using Fresa
using Profile
nx = 5
# specify the domain for the search, x ∈ [0,1]ⁿ
domain = Fresa.Domain(x -> zeros(length(x)), x -> ones(length(x)))
x = zeros(nx)
f = x -> ([ sum((x.-0.5).^2 .+ 1)
            sum(cos.(x))
            sum(sin.(x))],
          0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x,f)]
# now invoke Fresa to solve the problem
@profile for i=1
    pareto, population = Fresa.solve(f, p0, domain;
                                     archiveelite = false,
                                     npop=20, ngen=300,
                                     #output=100,
                                     tolerance=0.01)

    println("*** Pareto front:")
    println(population[pareto])
end
println("*** profile data")
println(": this may take some time so please wait")
Profile.print(format=:flat, sortedby=:count)
# testmultiobjective3 ends here
