# [[file:../fresa.org::testrosenbrock][testrosenbrock]]
using Fresa
# define the dimension of the problem and the initial guess
nx = 2
x0 = 0.5*ones(nx)
# specify the domain for the search, x ∈ [0,10]ⁿ
d = Fresa.Domain(x -> zeros(length(x)), x -> 10*ones(length(x)))
# the objective function is always feasible so needs only to return
# the actual objective function value
rosenbrock(x) = [(1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2]
# f = x -> ((x[1]-3)^2+(x[2]-5)^2+8, 0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x0,rosenbrock)]
# create a mutable structure that will be used to record information
# from the population analysis and create an instance of this
# structure
mutable struct BestFound
    gen
    point
end
global analysis
analysis = BestFound(0,missing)
# the function to analyse the population.  The arguments are described
# in the documentation for the Fresa.solve method.  The key thing to
# remember is that the solutions in the vectors nondom and dom are
# Fresa.Point objects.
function analyse(gen, nondom, nondomfit, dom, domfit)
    # if the current best is "better" than the previously saved best,
    # or a best solution has not yet been saved at all, save the new
    # solution as the best overall and record which generation this
    # happened.
    global analysis
    if analysis.point isa Missing || Fresa.dominates(nondom[1].z, analysis.point.z)
        analysis.gen = gen
        analysis.point = nondom[1]
    end
end
# the analysis function given to Fresa uses
# now invoke Fresa to solve the problem
best, pop = Fresa.solve(rosenbrock, p0;
                        analysePopulation = analyse,
                        domain=d,
                        ngen=1000,
                        issimilar = Fresa.similarx,
                        ϵ=1e-8)
println("Population at end: $pop")
println("Best solution is f($( best.x ))=$( best.z ) with g=$( best.g )")
println("which was found by generation $(analysis.gen)")
# testrosenbrock ends here
