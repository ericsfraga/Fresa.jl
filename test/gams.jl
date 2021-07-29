# [[file:../fresa.org::gamsfmo][gamsfmo]]
; function fmo(x::Array{Float64,1})
    open("gamsexample.gms", "w") do f
        write(f, "\$include gamsdeclarations.gms\n")
        write(f, "X1.fx = $(x[1]); \n")
        write(f, "X2.fx = $(x[2]); \n")
        write(f, "X3.fx = $(x[3]); \n")
        write(f, "solve TEST using NLP minimizing Z; \n")
        write(f, "file fresa /'gamsoutput.txt'/ ;\n")
        write(f, "put fresa ;\n")
        write(f, "put z.l /;\n")
        write(f, "put res.l /;\n")
        write(f, "put TEST.modelstat /;\n")
    end
    # execute GAMS
    run( `/opt/gams/latest/gams gamsexample.gms` )
    # read in results
    z = [0.0; 0.0]
    g = 0.0;
    open("gamsoutput.txt", "r") do f
        lines = readlines(f)
        z[1] = parse(Float64, lines[1])
        z[2] = abs(parse(Float64, lines[2]))
        modelstat = parse(Float64, lines[3])
        if modelstat != 1 && modelstat != 5
            g = 1
        end
    end
    # return results
    ( z, g )
end;
# gamsfmo ends here

# [[file:../fresa.org::gamsfsingle][gamsfsingle]]
; function fsingle(x::Array{Float64,1})
    open("gamsexample.gms", "w") do f
        write(f, "\$include gamsdeclarations.gms\n")
        write(f, "X1.fx = $(x[1]); \n")
        write(f, "X2.fx = $(x[2]); \n")
        write(f, "X3.fx = $(x[3]); \n")
        write(f, "solve TEST using NLP minimizing Z; \n")
        write(f, "file fresa /'gamsoutput.txt'/ ;\n")
        write(f, "put fresa ;\n")
        write(f, "put z.l /;\n")
        write(f, "put res.l /;\n")
        write(f, "put TEST.modelstat /;\n")
    end
    # execute GAMS
    run( `/opt/gams/latest/gams gamsexample.gms` )
    # read in results
    z = 0.0
    g = 0.0
    open("gamsoutput.txt", "r") do f
        lines = readlines(f)
        z = parse(Float64, lines[1])
        g = abs(parse(Float64, lines[2]))
        modelstat = parse(Float64, lines[3])
        if modelstat != 1 && modelstat != 5
            g = 10 # penalty function
        end
    end
    # return results
    ( z, g )
end;
# gamsfsingle ends here

# [[file:../fresa.org::*solve the multi-objective problem using Fresa][solve the multi-objective problem using Fresa:1]]
; using Fresa
domain = Fresa.Domain(x -> [0.0;0.0;0.0], x -> [5.0;3.0;3.0])
x0 = [4.0;2.0;2.0]
# create the initial population consisting of this single point
p0 = [Fresa.createpoint(x0,fmo)]
# now invoke Fresa to solve the problem
pareto, population = Fresa.solve(fmo, p0, domain;
                                 fitnesstype = :borda,
                                 ngen = 100)
println("Pareto front:")
println(population[pareto]);
# solve the multi-objective problem using Fresa:1 ends here

# [[file:../fresa.org::*solve the multi-objective problem using Fresa][solve the multi-objective problem using Fresa:2]]
; using PyPlot
z = [population[pareto[i]].z for i in 1:length(pareto)];
PyPlot.plot([z[i][1] for i=1:length(z)],
            [z[i][2] for i=1:length(z)],
            "ro")
PyPlot.savefig("gamsmo.pdf");
# solve the multi-objective problem using Fresa:2 ends here

# [[file:../fresa.org::*solve the single objective version][solve the single objective version:1]]
; best, pop = Fresa.solve(fsingle, p0, domain; ngen = 100)
println("Population: $pop")
println("Best: f($(best.x)) = $(best.z), $( best.g )");
# solve the single objective version:1 ends here
