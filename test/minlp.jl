# [[file:../fresa.org::testminlp][testminlp]]
using Distributed
using Printf
@everywhere using Fresa
# define new type for mixed integer problems
# in general, this would be vectors of real and integer values
@everywhere struct MI
    x :: Float64
    y :: Int32
end
import Base
Base.show(io::IO, m::MI) = print(io, m.x, " ", m.y)
f = s -> (3s.y - 5s.x,
          max(2s.y + 3s.x - 24,
              3s.x - 2s.y - 8,
              2s.y^2 - 2*√s.y + 11s.y + 8s.x - 39 - 2*√s.x*s.y^2))
# bounds
domain = Fresa.Domain(x -> MI(1.0, 1),
                      x -> MI(6.0, 6))
# function to find a neighbouring solution for MI type decision points
function Fresa.neighbour(s :: MI,
                         a :: MI,
                         b :: MI,
                         f :: Float64) :: MI
    x = s.x + (b.x-a.x)*(1-f)*2*(rand()-0.5)
    x = x < a.x ? a.x : (x > b.x ? b.x : x)
    # for the integer variable, we move in one direction or the other
    # a random number of places depending on fitness
    positive = rand(Bool)
    r = rand()
    # @printf(": neighbour: f=%g r=%g\n", f, r)
    inc = ceil(f*r*(b.y-a.y)/2)
    # @printf(": neighbour: positive=%s inc=%d\n", positive, inc)
    y = s.y + (positive ? inc : -inc)
    y = y < a.y ? a.y : (y > b.y ? b.y : y)
    return MI(x,y)
end
# create the initial population consisting of a single MI point
p0 = [Fresa.createpoint(MI(1.0, 1),f)]
# now invoke Fresa to solve the problem
best, pop = Fresa.solve(f, p0, domain; ngen=100)
println("Population: $pop")
println("Best: f($(best.x)) = $(best.z), $(best.g)")
# testminlp ends here

# [[file:../fresa.org::testminlpsupplement][testminlpsupplement]]
println("#+plot: ind:3 deps:(2) with:\"linespoints pt 7\" set:nokey set:\"yrange [0:1]\"")
ancestor = best.ancestor;
while ancestor != Some(nothing) && typeof(ancestor) != Nothing
    global ancestor
    println("| $(ancestor.point.z) | $(ancestor.fitness) | $(ancestor.generation) |")
    ancestor = ancestor.point.ancestor
end
# testminlpsupplement ends here
