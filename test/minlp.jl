# [[file:../fresa.org::*Example 1: Westerlund & Westerlund][Example 1: Westerlund & Westerlund:3]]
using Fresa
struct MI
    x :: Vector{Float64}
    y :: Vector{Int}
end
    function Fresa.neighbour(s :: MI,
                             a :: MI,
                             b :: MI,
                             f :: Float64) :: MI
        # use the neighbour function in Fresa to find a neighbour for the
        # floating point numbers in the representation of the current
        # point
        x = Fresa.neighbour(s.x, a.x, b.x, f)
        # the integer variables are treated differently.  We only consider
        # changing any value at all if the random number is greater than
        # the fitness value, which means that the most fit solutions will
        # likely not have the integer values changed.  If one is to be
        # changed, we limit to just one variable at a time.
        y = copy(s.y)
        if rand() > f
            i = rand(1:length(y))
            # consider the case of binary variables as special cases:
            # toggle the boolean value (which is essentially what a binary
            # variable can be considered to be); otherwise, change value
            # up or down randomly.
            if a.y[i] == 0 && b.y[i] == 1
                # binary variable
                y[i] = 1 - y[i]
            else
                # for the integer variable we select to change, we move in one
                # direction or the other a random number of places depending
                # on fitness
                positive = rand(Bool)
                # random number to decide how much change to make to this
                # integer
                r = rand()
                inc = ceil(f*r*(y[i]-y[i])/2)
                # @printf(": neighbour: positive=%s inc=%d\n", positive, inc)
                y[i] = y[i] + (positive ? inc : -inc)
                # keep within bounds
                y[i] = y[i] < a.y[i] ? a.y[i] : (y[i] > b.y[i] ? b.y[i] : y[i])
            end
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
