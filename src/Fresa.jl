# [[file:../fresa.org::modulestart][modulestart]]
# All code copyright © Eric S Fraga. 
# Licence for use and sharing can be found at
#   https://github.com/ericsfraga/Fresa.jl/blob/master/LICENSE
# Date of last change in version variable below.
module Fresa
# modulestart ends here

# [[file:../fresa.org::init][init]]
version = "8.0.1"
lastchange = "[2023-06-02 12:30+0100]"
using Dates                     # for org mode dates
using LinearAlgebra             # for norm function
using Printf                    # for formatted output
function __init__()
    println("# Fresa 🍓 PPA v$version, last change $lastchange")
end
# init ends here

# [[file:../fresa.org::pointtype][pointtype]]
"""

Point (`x`) in the search space along with objective function values
(`z[]`) and feasbility indication (`g`).  The type of `x` is problem
specific.  `z[]` and `g` hold `Float64` values.  `g` should be of
length 1.

"""
struct Point
    x :: Any                    # decision point
    z :: Vector                 # objective function values
    g :: Float64                # constraint violation
    ancestor                    # the parent of this point
    function Point(x,           # point in search space
                   f,           # objective function 
                   parameters = nothing, # arguments to objective function 
                   ancestor = nothing)   # for analysis of search process
        z = 0
        g = 0
        if ! ( parameters isa Nothing )
            (z, g) = f(x, parameters)
        else
            (z, g) = f(x)
        end
        if g isa Int
            g = float(g)
        end
        p = Nothing
        if rank(z) == 1
            p = new(x, z, g, ancestor)
        elseif rank(z) == 0
            p = new(x, [z], g, ancestor)
        else
            error("Fresa can only handle scalar and vector criteria, not $(typeof(z)).")
        end
        return p
    end
end
# pointtype ends here

# [[file:../fresa.org::showpoint][showpoint]]
import Base
Base.show(io::IO, p::Fresa.Point) = print(io, "f(", p.x, ")=", p.z, " g=", p.g)
# and also an array of points
function Base.show(io::IO, p::Array{Point,1})
    np = length(p)
    if np > 0
        nz = length(p[1].z)
        println(io, "|-")
        for i=1:nz
            print(io,"| z$(i) ")
        end
        println(io, "| g | x |")
        println(io,"|-")
        for i=1:length(p)
            for j=1:nz
                print(io,"| ", p[i].z[j], " ")
            end
            print(io, "| ", p[i].g, " ")
            print(io, "| ", p[i].x, " |\n")
        end
        println(io,"|-")
    else
        print(io,"empty")
    end
end
# showpoint ends here

# [[file:../fresa.org::pointsize][pointsize]]
import Base.size
Base.size(p :: Point) = ()
# pointsize ends here

# [[file:../fresa.org::ancestortype][ancestortype]]
struct Ancestor
    point :: Point        # the actual ancestor point
    fitness :: Float64    # the fitness of the ancestor
    generation :: Int32   # the generation when this point was created
end
# ancestortype ends here

# [[file:../fresa.org::*Ancestor <<ancestor>>][Ancestor <<ancestor>>:2]]
ancestor(p :: Point) = p.ancestor :: Union{Ancestor,Nothing}
# Ancestor <<ancestor>>:2 ends here

# [[file:../fresa.org::domaintype][domaintype]]
struct Domain
    lower                       # function which returns lower bound on search variable(s)
    upper                       # function which returns upper bound on search variable(s)
end
# domaintype ends here

# [[file:../fresa.org::createpoint][createpoint]]
function createpoint(x,f,parameters = nothing,ancestor = nothing)
    z = 0
    g = 0
    if ! ( parameters isa Nothing )
        (z, g) = f(x, parameters)
    else
        (z, g) = f(x)
    end
    if g isa Int
        g = float(g)
    end
    p = Nothing
    if rank(z) == 1
        p = Point(x, z, g, ancestor)
    elseif rank(z) == 0
        p = Point(x, [z], g, ancestor)
    else
        error("Fresa can only handle scalar and vector criteria, not $(typeof(z)).")
    end
    return p
end
# createpoint ends here

# [[file:../fresa.org::fitness][fitness]]
function fitness(pop, fitnesstype, steepness, generation, ngen)
    l = length(pop)
    indexfeasible = (1:l)[map(p->p.g,pop) .<= 0]
    indexinfeasible = (1:l)[map(p->p.g,pop) .> 0]
    @debug "Feasible/infeasible breakdown" indexfeasible indexinfeasible maxlog=3
    fit = zeros(l)
    factor = 1              # for placement in fitness interval (0,1)
    if length(indexfeasible) > 0
        feasible = view(pop,indexfeasible)
        # use objective function value(s) for ranking
        feasiblefit = vectorfitness(map(p->p.z,feasible), fitnesstype, steepness, generation, ngen)
        if length(indexinfeasible) > 0
            feasiblefit = feasiblefit./2 .+ 0.5 # upper half of fitness interval
            factor = 2                        # have both feasible & infeasible
        end
        fit[indexfeasible] = (feasiblefit.+factor.-1)./factor
    end
    if length(indexinfeasible) > 0
        # squeeze infeasible fitness values into (0,0.5) or (0,1) depending
        # on factor, i.e. whether there are any feasible solutions as well or not
        infeasible = view(pop,indexinfeasible)
        # use constraint violation for ranking as objective function values
        # may not make any sense given that points are infeasible
        fit[indexinfeasible] = vectorfitness(map(p->p.g, infeasible),
                                             fitnesstype,
                                             steepness,
                                             generation,
                                             ngen
                                             ) / factor;
    end
    fit
end
# fitness ends here

# [[file:../fresa.org::vectorfitness][vectorfitness]]
"""
For single objective problems, the fitness is simply the normalised
objective function value.

For multi-objective cases, there are three alternative measures of
fitness ranking possible.  The first is based on the Hadamard product
of the rank of each member of population accoring to each
criterion.  The second is based on a weighted Borda ranking based on
each criterion ranking.  Finally, a measure based on dominance,
similar to that used by the popular NSGA-II genetic algorithm, is
available.

"""
function vectorfitness(v, fitnesstype, steepness, generation, ngen)
    # determine number of objectives (or pseudo-objectives) to consider in
    # ranking
    l = length(v)
    if l == 1
        # no point in doing much as there is only one solution
        fit = [0.5]
    else
        m = length(v[1])
        # println("VF: v=$v")
        # println("  : of size $(size(v))")
        if m == 1                   # single objective 
            fitness = [v[i][1] for i=1:l]
        else                  # multi-objective
            rank = ones(m,l); #rank of each solution for each objective function 
            if fitnesstype == :hadamard
                for i=1:m
                    rank[i,sortperm([v[j][i] for j=1:l])] = 1:l
                end
                # hadamard product of ranks
                fitness = map(x->prod(x), rank[:,i] for i=1:l)
            elseif fitnesstype == :borda
                for i=1:m
                    rank[i,sortperm([v[j][i] for j=1:l])] = 1:l
                end
                # borda sum of ranks
                fitness = map(x->sum(x), rank[:,i] for i=1:l)
            elseif fitnesstype == :nondominated
                # similar to that used by NSGA-II (Deb 2000)
                fitness = zeros(l)
                maxl = assigndominancefitness!(fitness,v,1)
                # println("Resulting fitness: $fitness")
            else
                throw(ArgumentError("Type of fitness evaluation must be either :borda, :nondominated, or :hadamard, not $(repr(fitnesstype))."))
            end
        end
        # normalise (1=best, 0=worst) while avoiding
        # extreme 0,1 values using the hyperbolic tangent
        fit = adjustfitness(fitness, steepness, generation, ngen)
        # println(":  scaled fitness: $fit")
        @debug "Fitness calculations" v[1][1] v[2][1] v[l][1] fitness[1] fitness[2] fitness[l] fit[1] fit[2] fit[l] maxlog=3
    end
    fit
end
# vectorfitness ends here

# [[file:../fresa.org::adjustfitness][adjustfitness]]
function adjustfitness(fitness, steepness, generation, ngen)
    if (maximum(fitness)-minimum(fitness)) > eps()
        s = steepness
        if steepness isa Tuple
            a = (2*steepness[1]-2*steepness[2])/3
            b = - (3*steepness[1] - 3*steepness[2])/ngen^2
            d = steepness[1]
            s = a*generation^3 + b*generation^2 + c*generation + d
            @debug "Steepness " s "at generation" g
        end  
        fit = 0.5*(tanh.(4*s*(maximum(fitness) .- fitness)
                         / (maximum(fitness)-minimum(fitness))
                         .- 2*s) .+ 1)
    else
        # only one solution (or all solutions the same) in population
        fit = 0.5*ones(length(fitness))
    end
    fit
end
# adjustfitness ends here

# [[file:../fresa.org::assigndominancefitness][assigndominancefitness]]
function assigndominancefitness!(f,v,l)
    # assign value l to all members of v which dominate rest and then
    # recurse on those which are dominated
    (p, d) = paretoindices(v)
    # println("Assigning fitness $l to $p")
    f[p] .= l
    if !isempty(d)
        assigndominancefitness!(view(f,d),v[d],l+1)
    else
        l
    end
end
# assigndominancefitness ends here

# [[file:../fresa.org::neighbourarray][neighbourarray]]
function neighbour(x :: Vector{Float64},
                   f :: Float64,
                   d :: Domain
                   ) :: Vector{Float64}
    # allow movements both up and down in the domain for this variable
    # so determine the actual domain lower and upper bounds
    a = d.lower(x)
    b = d.upper(x)
    xnew = x .+ (1.0 .- f) .* 2(rand(length(x)).-0.5) .* (b.-a)
    xnew[xnew.<a] = a[xnew.<a];
    xnew[xnew.>b] = b[xnew.>b];
    return xnew
end
# neighbourarray ends here

# [[file:../fresa.org::neighbourfloat][neighbourfloat]]
function neighbour(x :: Float64,
                   f :: Float64,
                   d :: Domain
                   ) :: Float64
    # allow movements both up and down
    # in the domain for this variable
    a = d.lower(x)
    b = d.upper(x)
    newx = x + (b-a)*(2*rand()-1)/2.0 * (1-f)
    if newx < a
        newx = a
    elseif newx > b
        newx = b
    end
    newx
end
# neighbourfloat ends here

# [[file:../fresa.org::dominates][dominates]]
function dominates(a, b)
    all(a .<= b) && any(a .< b)
end
≻(a,b) = dominates(a,b)
# dominates ends here

# [[file:../fresa.org::*find Pareto set][find Pareto set:1]]
function paretoindices(z)
    n = length(z)
    dominance = [reduce(&, [!(z[i] ≻ z[j]) for i ∈ 1:n]) for j ∈ 1:n]
    paretoindices = filter(j -> dominance[j], 1:n)
    dominatedindices = filter(j -> !dominance[j], 1:n)
    (paretoindices, dominatedindices)
end
# find Pareto set:1 ends here

# [[file:../fresa.org::pareto][pareto]]
# indices of non-dominated and dominated points from the population of
# Point objects
function pareto(pop :: Vector{Point})
    l = length(pop)
    indexfeasible = (1:l)[map(p->p.g,pop) .<= 0]
    indexinfeasible = (1:l)[map(p->p.g,pop) .> 0]
    if length(indexfeasible) > 0
        subset = view(pop,indexfeasible)
        indices = indexfeasible
    else
        #println(": Fresa.pareto warning: no feasible solutions.  Pareto set meaningless?")
        subset = pop
        indices = 1:l
    end
    z = map(p->p.z, subset)
    # use function below to return indices of non-dominated and
    # dominated from objective function values alone in the subset of
    # feasible solutions
    (p, d) = paretoindices(z)
    (indices[p], indices[d])
end
# pareto ends here

# [[file:../fresa.org::printhistorytrace][printhistorytrace]]
function printHistoryTrace(p :: Point)
    a = p.ancestor
    while ! (a isa Nothing)
        println("| $(a.generation) | $(a.fitness) |")
        a = a.point.ancestor
    end
end
# printhistorytrace ends here

# [[file:../fresa.org::prune][prune]]
function prune(pop :: AbstractArray, issimilar, ϵ, domain)
    l = length(pop)
    # we will return a diverse population where similar solutions have
    # been removed
    diverse = [pop[1]]
    # consider each solution in the population
    for i=2:l
        similar = false
        j = 0
        # compare this solution with all already identified as diverse
        # enough
        while !similar && j < length(diverse)
            j += 1
            similar = issimilar(diverse[j], pop[i], ϵ, domain)
        end
        if !similar
            push!(diverse,pop[i])
        end
    end
    # return diverse population and count of points removed
    (diverse, length(pop)-length(diverse))
end
# prune ends here

# [[file:../fresa.org::similarx][similarx]]
function similarx(p1, p2, ϵ, domain)
    norm(p1.x-p2.x) < ϵ &&      # decision variables
        norm(p1.g-p2.g) < ϵ &&  # difference in violation
        ( (p1.g ≤ 0 && p2.g ≤ 0) || (p1.g > 0 && p2.g > 0)) # both same feasibility
end
# similarx ends here

# [[file:../fresa.org::similarz][similarz]]
function similarz(p1, p2, ϵ, domain)
    norm(p1.z-p2.z) < ϵ &&      # objective function values
        norm(p1.g-p2.g) < ϵ &&  # difference in violation
        ( (p1.g ≤ 0 && p2.g ≤ 0) || (p1.g > 0 && p2.g > 0)) # both same feasibility
end
# similarz ends here

# [[file:../fresa.org::randompopulation][randompopulation]]
function randompopulation(n, f, parameters, p0, domain :: Domain)
    p = Point[]                 # population object
    for j in 1:n
        # l = domain.lower(p0.x)
        # @show l
        # u = domain.upper(p0.x)
        # @show u
        # x = randompoint(l,u)
        # push!(p, createpoint(x, f, parameters))
        push!(p, Point(randompoint(domain.lower(p0.x), domain.upper(p0.x)),
                       f, parameters))
    end
    p
end
# randompopulation ends here

# [[file:../fresa.org::randompoint][randompoint]]
function randompoint(a :: Float64, b :: Float64)
    x = a + rand()*(b-a)
end
function randompoint(a, b)
    x = a + rand(length(a)).*(b-a)
end
# randompoint ends here

# [[file:../fresa.org::select][select]]
function select(f, size)
    indices = rand(1:length(f), size)       # generate size indices
    best = argmax([f[i] for i ∈ indices])
    indices[best]
end
# select ends here

# [[file:../fresa.org::solve][solve]]
""" 

Solve an optimisation problem, defined as the minimization of the
values returned by the objective function, `f`.  `f` returns not only
the objective function values, an array of `Float64` values, but also
a measure of feasibility (≤0) or infeasibility (>0).  The problem is
solved using the Fresa algorithm.  `p0` is the initial population
which has to have at least one member, a `Point`, and `domain`
describes the search domain.  This latter argument is an instance of
the `Fresa.Domain` struct which has a `lower` and an `upper` members
which are functions to be evaluated with a current point in the
domain.

The return values for the solution of a single criterion problem are
the best point and the full population at the end of the search. 

For a multi-objective problem, the returned values are the set of
indices for the points within the full population (the second returned
value) approximating the *Pareto* front.

The population will consist of an array of `Fresa.Point` objects, each
of which will have the point in the search space, the objective
function value and the feasibility measure.

"""
function solve(f, p0;                # required arguments
               archiveelite = false, # save thinned out elite members
               domain = nothing,     # search domain: will often be required but not always
               elite = true,         # elitism by default
               ϵ = 0.0001,           # ϵ for similarity detection
               fitnesstype = :hadamard, # how to rank solutions in multi-objective case
               issimilar = nothing,  # function for diversity check: see prune function
               multithreading = false, # use multiple threads for objective function evaluation
               ngen = 0,             # stopping criterion: number of generations
               np = 10,              # points to propagate: constant (single value) or dynamic (tuple)
               nfmax = 0,            # stopping criterion: maximum number of function evaluations
               nrmax = 5,            # number of runners maximum
               ns = 100,             # number of stable solutions for stopping
               orglevel = "",        # default org mode heading indentation for any output
               output = 1,           # how often to output information
               parameters = nothing, # allow parameters for objective function 
               plotvectors = false,  # generate output file for search plot
               populationoutput = false, # output population every generation?
               tournamentsize = 2,   # number to base selection on
               steepness = 1.0,      # show steep is the adjustment shape for fitness
               usemultiproc = false) # parallel processing by Fresa itself?
    output > 0 && println("$orglevel* Fresa solve $f $(orgtimestamp(now()))")
    tstart = time()
    nf = 1                   # number of function evaluations
    npruned = 0              # number solutions pruned from population
    nz = length(p0[1].z)     # number of criteria
    pop = copy(p0);          # create/initialise the population object
    if archiveelite
        archive = Point[]
    end
    # check to see if at least one stopping criterion has been
    # defined.  If not, set the number of generations to an arbitrary
    # value but notify the user
    if ngen ≤ 0 && nfmax ≤ 0
        @warn "No stopping criteria given so defaulting to ngen=100"
        ngen = 100
        nfmax = typemax(Int32)  # arbitrarily large?
    end
    if output > 0
        println("#+name: $(f)settings")
        println("| variable | value |")
        println("|-")
        println("| elite | $elite |")
        println("| archive | $archiveelite |")
        println("| ϵ | $ϵ |")
        println("| fitness | $fitnesstype |")
        println("| issimilar | $issimilar |")
        println("| multithreading | $multithreading |")
        if nfmax ≤ 0
            println("| nfmax | ∞ |")
        else
            println("| nfmax | $nfmax |")
        end
        if ngen ≤ 0
            println("| ngen | ∞ |")
        else
            println("| ngen | $ngen |")
        end
        println("| np | $np |")
        println("| nrmax | $nrmax |")
        println("| ns | $ns |")
        println("| steepness | $steepness |")
        println("| tournamentsize | $tournamentsize |")
        println("|-")
        # output != 0 && println(": solving with ngen=$ngen np=$np nrmax=$nrmax ns=$ns")
        # output != 0 && println(": elite=$elite archive elite=$archiveelite fitness type=$fitnesstype")
    end
    if plotvectors
        plotvectorio = open("fresa-vectors-$(orgtimestamp(now())).data", create=true, write=true)
        output > 0 && println(": output of vectors for subsequent plotting")
    end
    # if np was given as a tuple, we are to have a dynamic
    # population size.  This only makes sense for multi-objective
    # optimization problems so a warning will be given otherwise.
    npmin = np
    npmax = np
    if np isa Tuple
        if nz > 1
            npmin = np[1]
            npmax = np[2]
            if npmin > npmax
                error("Dynamic population sizing requires min <= max; you specified $np")
            end
            np = npmin      # start with minimum possible
        else
            println("*Warning*: you have specified a tuple for population size: $np")
            println("This only makes sense for multi-objective optimization problems.")
            println("np will be set to $(np[1]).")
            np = np[1]      # be optimistic and use minimum given
        end
    end
    # we use multithreading if asked for *and* if we have more than
    # one thread available
    multithreading = multithreading && Threads.nthreads() > 1 
    # we use parallel computing if we have more than one processor
    parallel = usemultiproc && nprocs() > 1
    # parallel = false
    if output > 0
        println(": function evaluations performed ",
                parallel
                ? "in parallel with $(nprocs()) processors."
                : (multithreading
                   ? "in parallel with $(Threads.nthreads()) threads."
                   : "sequentially."))
        println("$orglevel** initial population")
        println("#+name: $(f)initial")
        println(pop)
    end
    if output > 0
        println("$orglevel** evolution")
        println("#+name: $(f)evolution")
        println("#+plot: ind:1 deps:(6) with:\"points pt 7\" set:\"logscale x\"")
        @printf("| %9s | %9s | %9s | %9s | %9s |", "gen", "pop",
                (elite && nz > 1) ? "pareto" : "nf", "pruned", "t (s)")
        for i in 1:nz
            @printf(" z%-8d |", i)
        end
        @printf(" %9s |", "g")
        @printf("\n|-\n")
    end
    # now evolve the population.  There are two potential stopping
    # criteria: a maximum number of generations to perform (ngen) and
    # a maximum number of function evaluations (nfmax).  Either one or
    # the other will determine how much work we do, whichever is
    # reached first.
    gen = 0
    while (ngen ≤ 0 || gen < ngen) && (nfmax ≤ 0 || nf < nfmax)
        gen += 1                # keep count of generations performed

        # evaluate fitness which is adjusted depending on value of
        # steepness, a value that may depend on the generation
        fit = fitness(pop, fitnesstype, steepness, gen, ngen)
        if gen == 1
            @debug "Initial fitness" f=fit
        end
        # sort
        index = sortperm(fit)
        if populationoutput
            println("\nGeneration $gen full population is:")
            println(pop)
            println("Fitness vector: $fit")
        end
        # and remember best which really only makes sense in single
        # criterion problems but is best in multi-objective case in
        # the ranking measure used by Fresa
        best = pop[index[end]]
        # if elitism is used
        if elite
            if nz > 1
                # elite set is whole pareto set unless it is too
                # big. Recall that the pareto function returns the set
                # of indices into the population
                wholepareto = pareto(pop)[1]
                # if using dynamic population sizing, adjust the population
                np = 2 * length(wholepareto)
                if np < npmin
                    np = npmin
                end
                if np > npmax
                    np = npmax
                end
                # now check that the pareto is not too big.  if it is, thin it out
                if length(wholepareto) > ceil(np/2)
                    newpop, removed = thinout(pop, fit, wholepareto, ceil(Int,np/2))
                    if archiveelite
                        # add removed solutions to the archive, pruning if desired
                        if issimilar != nothing
                            archive = prune(append!(archive, removed), issimilar, ϵ, domain)[1]
                        else
                            archive = append!(archive, removed)
                        end
                        # reduce archive to non-dominated solutions alone
                        archive = archive[pareto(archive)[1]]
                    end
                else
                    newpop = pop[wholepareto]
                end
            else
                # elite set is single element only
                newpop = [best]
            end
            # if plotting vectors for the search, include elitism
            if plotvectors
                for p in newpop
                    write(plotvectorio, "$(gen-1) $(p.x)\n$gen $(p.x)\n\n")
                end
            end
        else
            newpop = Point[]
        end
        if output >= 0
            print(stderr, ": $nf@$gen npop $(length(newpop))/$(length(pop))",
                  archiveelite ? " na=$(length(archive))" : "",
                  " most fit ",
                  best.g ≤ 0 ? "z=$(best.z)" : "g=$(best.g)",
                  " \r")
            # if output has been requested, check to see if output is
            # required now and then also check to see if the frequency
            # needs to be reduced.
            if output > 0
                if gen%output == 0 || gen == ngen
                    @printf("| %9d | %9d | %9d | %9d | %9.2f |", gen, length(fit),
                            (elite && nz > 1) ? length(newpop) : nf, npruned, time()-tstart)
                    for i = 1:length(best.z)
                        print(" $(best.z[i]) |")
                    end
                    print(" $(best.g) |")
                    println()
                end
                if 10^(floor(log10(gen))) > output
                    output = 10^(Int(floor(log10(gen))))
                end
            end
        end
        # if we are using any form of multiprocessing, either threads
        # or multiple cores, create an array to store all new points
        # which we evaluate later in parallel.  Ideally, also keep
        # track of the points from which new points are derived to
        # provide the backward link through the evolution but this is
        # currently disabled as the creation of the Ancestor object
        # requires more information than I am currently storing away.
        if multithreading || parallel
            x = Any[] # typeof(newpop[1].x)[]
            # points = Point[]
        end
        # now loop through population, applying selection and then
        # generating neighbours
        l = length(pop)
        for i in 1:min(l,np)
            s = select(fit, tournamentsize)
            # println(": selection $i is $s")
            # println(": size of pop is $(size(pop))")
            selected = pop[s]
            if !elite
                # if no elitism, we ensure selected members remain in population
                push!(newpop, selected)
                if plotvectors
                    write(plotvectorio, "$(gen-1) $(selected.x)\n$gen $(selected.x)\n\n")
                end
            end
            # number of runners to generate, function of fitness
            nr = ceil(fit[s]*nrmax*rand())
            if nr < 1
                nr = 1
            end
            # println(": generating $nr runners")
            for r in 1:nr
                # create a neighbour, also function of fitness,
                # optionally passing a Domain object for the search
                # space.
                if domain isa Nothing
                    newx = neighbour(pop[s].x, fit[s])
                else
                    newx = neighbour(pop[s].x, fit[s], domain)
                end
                nf += 1
                # for parallel evaluation, we store the neighbours and
                # evaluate them later; otherwise, we evaluate
                # immediately and save the resulting point
                if multithreading || parallel
                    push!(x, newx)
                    # push!(points, pop[s])
                else
                    push!(newpop, Point(newx, f, parameters, Ancestor(pop[s],fit[s],gen)))
                    if plotvectors
                        write(plotvectorio, "$(gen-1) $(pop[s].x)\n$gen $newx\n\n")
                    end
                end
            end
            # remove selected member from the original population so
            # it is not selected again
            splice!(fit, s)
            splice!(pop, s)
        end
        # if we are making use of parallel computing, we evaluate all
        # points generated in previous loop.  Parallel processing is
        # done either via multithreading or with multiple
        # processors.  The former is easier as it's based on shared
        # memory.
        if multithreading       # using threads and shared memory
            results = Array{Point}(undef,length(x))
            Threads.@threads for i ∈ 1:length(x)
                results[i] = Point(x[i],f,parameters)
            end
            append!(newpop, results)
            # elseif parallel        # using multiple processors with remote calls
            #     # will be used to collect results from worker processors
            #     results = Array{Future,1}(undef, nprocs())
            #     i = 0;
            #     while i < length(x)
            #         # issue remote evaluation call
            #         for j=1:nprocs()
            #             if i+j <= length(x) 
            #                 # TODO: the information about the ancestor is
            #                 # not available; this needs to be stored above
            #                 results[j] = @spawn createpoint(x[i+j],f,parameters)
            #                 nf += 1
            #             end
            #         end
            #         # now wait for results
            #         for j=1:nprocs()
            #             if i+j <= length(x)
            #                 push!(newpop, fetch(results[j]))
            #             end
            #         end
            #         i += nprocs()
            #     end
        end
        # and finally, if diversity control has been enabled
        # (issimilar function provided and tolerance specified),
        # remove any duplicate points in the new population and make
        # it the current population for the next generation;
        # otherwise, simply copy over
        if issimilar != nothing && ϵ > eps()
            (pop, nn) = prune(newpop, issimilar, ϵ, domain)
            npruned += nn
        else
            pop = newpop
        end
    end
    # if output has been requested, check to see if output is
    # required for the last generation
    if output > 0
        if gen%output != 0
            @printf("| %9d | %9d | %9d | %9d | %9.2f |", gen, length(fit),
                    (elite && nz > 1) ? length(pop) : nf, npruned, time()-tstart)
            println()
        end
    end
    output > 0 && println("$orglevel** Fresa run finished\n: nf=$nf npruned=$npruned", archiveelite ? " archived=$(length(archive))" : "")
    if plotvectors
        close(plotvectorio)
    end
    if nz == 1
        fit = fitness(pop, fitnesstype, steepness, ngen, ngen)
        index = sortperm(fit)
        best = pop[index[end]]
        return best, pop
    else
        return pareto(archiveelite ? append!(pop,archive) : pop)[1], pop
    end
end
# solve ends here

# [[file:../fresa.org::thinout][thinout]]
function thinout(pop, fit, pareto, n::Int)
    indices = sortperm(fit[pareto])
    return pop[pareto[indices[end-n+1:end]]], pop[pareto[indices[1:end-n]]]
end
# thinout ends here

# [[file:../fresa.org::orgtimestamp][orgtimestamp]]
function orgtimestamp(dt::DateTime)
    return @sprintf("[%d-%02d-%02d %02d:%02d]",
                    Dates.year(dt),
                    Dates.month(dt),
                    Dates.day(dt),
                    Dates.hour(dt),
                    Dates.minute(dt))
end
# orgtimestamp ends here

# [[file:../fresa.org::rank][rank]]
rank(x :: Any) = length(size(x))
# rank ends here

# [[file:../fresa.org::moduleend][moduleend]]
end
# moduleend ends here
