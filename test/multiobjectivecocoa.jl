# [[file:../fresa.org::testmultiobjectivecocoa][testmultiobjectivecocoa]]
using Cocoa
using Fresa
nx = 2
# specify the domain for the search, x ∈ [0,1]ⁿ
lower = zeros(nx)
upper = ones(nx)
# initial point in domain
x = rand(nx)
# objective function 
f = x -> ( [sin(x[1]-x[2]); cos(x[1]+x[2])], 0)
# create the initial population consisting of this single point
p0 = [Fresa.Point(x,f)]
# define the solvers that Cocoa will use.  In this case, we are going
# to consider two instances of Fresa, each with a different fitness
# calculation method.
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
solvers = [Cocoa.Solver("hadamard", hadamard)
           Cocoa.Solver("nondominated", nondominated)]
model = Cocoa.Model("mo", f, p0, lower, upper, nothing)
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
    println(best)
end
# testmultiobjectivecocoa ends here
