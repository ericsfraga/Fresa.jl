# [[file:../fresa.org::testmultiobjective3][testmultiobjective3]]
using Fresa
using Profile
nx = 5
a = zeros(nx)
b = ones(nx)
x = zeros(nx)
f = x -> ([ sum((x.-0.5).^2 .+ 1)
            sum(cos.(x))
            sum(sin.(x))],
          0)
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x,f)]
# now invoke Fresa to solve the problem
@profile for i=1
    pareto, population = Fresa.solve(f, p0, a, b;
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
