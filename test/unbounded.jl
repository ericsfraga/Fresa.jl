# [[file:../fresa.org::unboundedtest][unboundedtest]]
# first thing to do is create a test objective function.  This
# function will not be challenging in a traditional sense, as it is
# simply a two dimensional quadratic, but the minimum is shifted by an
# arbitrary distance in both dimensions meaning that the search domain
# is not known a priori.
α = tan((π - eps()) * (rand()-0.5))
β = tan((π - eps()) * (rand()-0.5))
println("Using α=$α and β=$β")
f(x) = ( (x[1]-α)^2 + (x[2]-β)^2, 0.0 )
# now try to solve this
x0 = [0.0, 0.0]
f(x0)
using Fresa
p0 = [Fresa.Point(x0,f)]
best, pop = Fresa.solve(f, p0)
println("Best obtained: $best")
println(pop)
# unboundedtest ends here
