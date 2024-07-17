# [[file:../fresa.org::testmultiobjectivecocoa][testmultiobjectivecocoa]]
using Cocoa
using Fresa
# the size of the problem: number of real decision variables
nx = 2
# specify the domain for the search, x ∈ [0,1]ⁿ
lower = zeros(nx)
upper = ones(nx)
# random initial point in domain
x = rand(nx)
# the bicriteria objective function which is always feasible 
f = x -> ( [sin(x[1]-x[2]); cos(x[1]+x[2])], 0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x,f)]

#=

,----
| Start of Cocoa specific configuration.
`----

If this were just Fresa, we would solve the problem by invoking the
Fresa.solve() method; instead, we define the solvers that Cocoa will
be able to use, then inform Cocoa about these solvers, define the
Model for Cocoa, and finally invoke the Cocoa.solve() method.

=#

# define the solvers that Cocoa will use.  In this case, we are going
# to consider two instances of Fresa, each with a different fitness
# calculation method.  The calling sequence for a solver from Cocoa is
# always just an initial population, which will have been passed to
# Cocoa itself, the search domain defined by lower and upper bounds,
# and some Cocoa specific parameters which *must* be passed on to the
# individual solver, Fresa.solve() in this case.

# the two solvers we define differ purely in the fitness method used
# to rank multi-objective solutions.

# note that the stopping criterion specified should be large enough to
# never be met.  Cocoa itself will have a stopping criterion that will
# terminate the search.  Also, Fresa should be told to not generate
# any output as Cocoa should output all relevant progress information
# during the search.

# finally, Fresa must be told that it is running from within Cocoa by
# setting the cocoa argument to true.

hadamard(p0, lower, upper, parameters) = Fresa.solve(Cocoa.objectivefunction,
                                                     p0,
                                                     cocoa = true,
                                                     fitnesstype = :hadamard,
                                                     lower = lower,
                                                     nfmax = 1e6,
                                                     output = 0,
                                                     parameters = parameters,
                                                     upper = upper)

nondominated(p0, lower, upper, parameters) = Fresa.solve(Cocoa.objectivefunction,
                                                         p0,
                                                         cocoa = true,
                                                         fitnesstype = :nondominated,
                                                         lower = lower,
                                                         nfmax = 1e6,
                                                         output = 0,
                                                         parameters = parameters,
                                                         upper = upper)

# once the solvers have been defined, we create a vector of
# Cocoa.Solver types where each solver is defined by a name and the
# function defined above.

solvers = [Cocoa.Solver("hadamard", hadamard)
           Cocoa.Solver("nondominated", nondominated)]

# the Model for Cocoa is the objective function, f in this case, the
# initial population, the search domain (lower and upper bounds) and
# any parameters (nothing in this example) that may be required for
# the evaluation of the objective function.

model = Cocoa.Model("mo", f, p0, lower, upper, nothing)

# we are now ready to invoke the Cocoa solve method.  However, Cocoa
# supports cooperation between solvers (and competition but that's
# left for another day).  The features argument to Cocoa.solve()
# allows us to indicate whether we want the solvers to cooperate by
# having any good solution found by one shared with other solvers.  We
# consider solving the problem once without sharing and once with
# sharing.

# define a sequence of tests for features in Cocoa
tests = [Set()
         Set([:sharebest])]
for test ∈ tests
    println("* test features $test")
    best = Cocoa.solve(model, solvers,
                       features = test,
                       nm = 2,
                       nt=1_000,
                       orglevel = "**",
                       output = true)
    # best will be a vector of Cocoa solution types, each of which is
    # non-dominated.  A Cocoa solution includes the design point, d,
    # the objective function values, z, and the infeasibility measure,
    # g.  We can create Fresa.Point types for each of these:
    population = [Fresa.Point(s.d, s.z, s.g) for s in best]
    println(population)
end
# testmultiobjectivecocoa ends here
