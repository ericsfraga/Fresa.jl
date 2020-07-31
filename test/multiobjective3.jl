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
@profile for i=1
    pareto, population = Fresa.solve(f, x, a, b;
                                     archiveelite = false,
                                     npop=20, ngen=300,
                                     #output=100,
                                     tolerance=0.01)
    println("*** profile data")
    Profile.print(format=:flat, sortedby=:count)

    println("*** Pareto front:")
    println(population[pareto])

    using PyPlot
    z = [population[pareto[i]].z for i in 1:length(pareto)];
    PyPlot.plot3D([z[i][1] for i=1:length(z)],
                  [z[i][2] for i=1:length(z)],
                  [z[i][3] for i=1:length(z)],
                  "ro")
    PyPlot.savefig("x3.pdf")
end
# testmultiobjective3 ends here
