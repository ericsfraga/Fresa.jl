# [[file:../fresa.org::*Example 2: Quesada & Grossmann][Example 2: Quesada & Grossmann:2]]
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
    f = s -> (10*s.x[1]^2 - s.x[2] + 5*(s.y[1] - 1),
              max(
                  s.x[2] - 5*log(s.x[1]+1) - 3*s.y[1],
                  s.x[1]^2 - s.x[2] - s.y[1] - 1,
                  s.x[1] + s.x[2] + 20*s.y[1] - 24,
                  3*s.x[1] + 2*s.x[2] - 10
              ))
    # bounds, avoiding the absolute lower bound at -1 on the x variables
    # as this will cause numeric difficulties with the log function. eps()
    # is a built-in function which is the machine precision.
    d = Fresa.Domain(x -> MI([-1.0+eps(), -1.0+eps()], [0]),
                     x -> MI([50.0, 50.0], [1]))
    # create the initial population consisting of a single MI point
    p0 = [Fresa.Point(MI([0.0, 0.0], [1]),f)]
    # now invoke Fresa to solve the problem
    best, pop = Fresa.solve(f, p0; domain=d, ngen=10_000, tolerance=0.0)
    println("Population: $pop")
    println("Best: f($(best.x)) = $(best.z), $(best.g)")
# Example 2: Quesada & Grossmann:2 ends here
