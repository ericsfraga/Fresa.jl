# [[file:../fresa.org::*Example 1: Westerlund & Westerlund][Example 1: Westerlund & Westerlund:3]]
using Fresa
struct MI
    x :: Vector{Float64}
    y :: Vector{Int}
end
    function Fresa.neighbour(s :: MI,
                             f :: Float64,
                             d :: Fresa.Domain) :: MI
        # find the lower and upper bounds on all variables
        a = d.lower(s)
        b = d.upper(s)
        # use the neighbour function in Fresa to find a neighbour for the
        # floating point numbers in the representation of the current
        # point; the domain is defined by the real parts of the overall
        # domain, a and b retrieved above.
        x = Fresa.neighbour(s.x, f, Fresa.Domain(x -> a.x, x -> b.x))
        # the integer variables are treated differently.  We only consider
        # changing any value at all if the random number is greater than
        # the fitness value, which means that the most fit solutions will
        # likely not have the integer values changed.  If one is to be
        # changed, we limit to just one variable at a time.
        y = rand() > f ?
            y = Fresa.neighbour(s.y, f, Fresa.Domain(s -> d.lower(s).y, s -> d.upper(s).y)) :
            copy(s.y)
        # return the new search point consisting of both the floating
        # point numbers and the integer numbers
        return MI(x,y)
    end
    # objective function and constraints
    f = s -> (3s.y[1] - 5s.x[1],
              max(2s.y[1] + 3s.x[1] - 24,
                  3s.x[1] - 2s.y[1] - 8,
                  2s.y[1]^2 - 2*√s.y[1] + 11s.y[1] + 8s.x[1] - 39 - 2*√s.x[1]*s.y[1]^2))
    # bounds
    d = Fresa.Domain(x -> MI([1.0], [1]),
                     x -> MI([6.0], [6]))
    # create the initial population consisting of a single MI point
    p0 = [Fresa.Point(MI([1.0], [1]),f)]
    # now invoke Fresa to solve the problem
    best, pop = Fresa.solve(f, p0; domain=d, ngen=100) # , populationoutput=true)
    println("Population: $pop")
    println("Best: f($(best.x)) = $(best.z), $(best.g)")
    println("#+plot: ind:3 deps:(2) with:\"linespoints pt 7\" set:nokey set:\"yrange [0:1]\"")
    ancestor = best.ancestor;
    while ancestor != Some(nothing) && ! (ancestor isa Nothing)
        global ancestor
        println("| $(ancestor.point.z) | $(ancestor.fitness) | $(ancestor.generation) |")
        ancestor = ancestor.point.ancestor
    end
# Example 1: Westerlund & Westerlund:3 ends here
