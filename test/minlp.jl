# [[file:../fresa.org::*Example 1: Westerlund & Westerlund][Example 1: Westerlund & Westerlund:3]]
using Fresa
struct MI
    x :: Vector{Float64}
    y :: Vector{Int32}
end
    function Fresa.neighbour(s :: MI,
                             a :: MI,
                             b :: MI,
                             f :: Float64) :: MI
        # use the neighbour function in Fresa to find a neighbour for the
        # floating point numbers in the representation of the current
        # point
        x = Fresa.neighbour(s.x, a.x, b.x, f)
        # now find a neighbour for the integer part: instead of possibly
        # changing all the integer values, as we do for the floating point
        # numbers, we consider changing a selection of the values with the
        # number to change based on the fitness: the better the fitness,
        # the fewer integer values we change.
        ni = ceil(length(s.y) * (1-f))
        y = copy(s.y)
        for i ∈ 1:ni
            j = rand(1:length(s.y))
            # for the integer variable we select to change, we move in one
            # direction or the other a random number of places depending
            # on fitness
            positive = rand(Bool)
            # random number to decide how much change to make to this
            # integer
            r = rand()
            inc = ceil(f*r*(b.y[j]-a.y[j])/2)
            # @printf(": neighbour: positive=%s inc=%d\n", positive, inc)
            y[j] = s.y[j] + (positive ? inc : -inc)
            # keep within bounds
            y[j]= y[j] < a.y[j] ? a.y[j] : (y[j] > b.y[j] ? b.y[j] : y[j])
        end
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
    domain = Fresa.Domain(x -> MI([1.0], [1]),
                          x -> MI([6.0], [6]))
    # create the initial population consisting of a single MI point
    p0 = [Fresa.createpoint(MI([1.0], [1]),f)]
    # now invoke Fresa to solve the problem
    best, pop = Fresa.solve(f, p0, domain; ngen=100)
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
