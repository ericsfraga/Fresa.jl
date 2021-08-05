---
title: 'Fresa: a plant propagation algorithm for black-box single and multiple objective optimization '
tags:
  - Julia
  - optimization
  - multi-objective
  - plant propagation algorithm
authors:
  - name: Eric S. Fraga
    orcid: 0000-0002-5732-6082
    affiliation: 1
affiliations:
 - name: Sargent Centre for Process Systems Engineerg, Department of Chemical Engineering, University College London (UCL)
   index: 1
date: 5 August 2021
bibliography: paper.bib

# Summary

`Fresa` implements a nature inspired plant propagation algorithm for the solution of single and multiple objective optimization problems. The method is population based and evolutionary.  Treating the objective function as a black box, the implementation is able to solve problems exhibiting behaviour that is challenging for mathematical programming methods.  Taking advantage of the dynamic typing and multiple dispatch capabilities of the Julia language `[@bezanson-etal-2017a]`, `Fresa` is easily adapted to new problems which may benefit from bespoke representations of solutions.  Further, the support for threads in Julia enables an efficient implementation.


# Statement of need

There are a number of different ways optimization methods can be categorized: mathematical programming versus direct search and deterministic versus stochastic being two examples of categories.  **Mathematical programming** methods typically require the model to be represented by a set of equations, equations that could be algebraic or differential, and make use of derivatives to guide a search to an optimum.  **Direct search** methods, otherwise known as *zeroth*-order methods, on the other hand, require only the objective function value itself, with possibly an indication of the feasibility of a given solution.  Mathematical programming methods are usually deterministic whereas direct search methods may be deterministic `[@kelley-1999]` or stochastic, such as genetic algorithms `[@goldberg-1989]`, simulated annealing `[@simulated-annealing-1987]`, and particle swarm optimization `[@488968]`, to mention only a few.  `Fresa` `[@fresa]` is an example of a stochastic direct search method.

Direct search methods, and especially those based on stochastic searches, are often necessary due to the properties of the optimization problem.  Problems with discontinuities, for instance, are challenging to handle robustly, if at all, using mathematical programming approaches.  Properties such as nonlinear and non-convex behaviour also challenge many solution methods as would models based on ordinary or partial differential equations.  In engineering, models may exhibit noise.  Noise can arise when the objective function is based on data from experiments but noise may also arise when the objective function requires the numerical solution of embedded differential equations, for instance.  Noise in the values of the objective function create difficulties in using derivative information for guiding the search and will, in the best case, lead to sub-optimal results.

Further, many problems, especially in design, have multiple criteria for the evaluation of points in the search space.  Multi-objective optimization is therefore a desirable feature of an optimization solver, although such problems can and are solved using single objective methods through multiple solutions of a suitably modified single objective problem (weighting the objectives; &epsilon;-constraint methods).


# Features

As noted above, `Fresa` is an example of a stochastic direct search method.  Specifically, it implements a population based evolutionary procedure.  New solutions in the search space are generated, using random numbers, with the identification of *neighbour* solutions inspired by the propagation of runners for Strawberry plants `[@salhi-fraga-2011a]`.  Solutions are selected for propagation based on their fitness which is related to the values of the objective function(s).  The number of runners and their distance for propagation are functions of the fitness as well.  `Fresa` supports both single and multi-objective optimization.

An advantage of the underlying algorithm is that it has few user-tunable parameters and the results have been shown to be relatively insensitive to the values of these parameters `[@dejonge-vandenberg-2020a]`.  These parameters include the population size, the number of cycles or generations to perform, the maximum number of runners to generate, and the choice of and parameters for the fitness method used.  

One interesting feature of `Fresa` is that it allows for heterogeneous populations where individuals may use different representations of potential solutions for the problem `[@fraga2021multiple]`.  This is achieved by the combination of dynamic typing and multiple dispatch provided by the Julia language.

For multi-objective optimization, different algorithms have been implemented to assign fitness values to individual solutions in the population.  One of these is the non-dominated sorting algorithm, exemplified by the popular NSGA-II genetic algorithm implementation `[@deb-2000]`.  In `Fresa`, other alternatives are fitness values based on the *Hadamard* product of the rankings with respect to the individual criteria `[@fraga-amusat-2016a]` and another alternative based on the *Borda* sum of those individual rankings.  These latter two fitness assignment methods tend to emphasise solutions found towards the ends of the Pareto frontier whereas the non-dominated sorting algorithm may lead to the ends of the frontier having less representation in the population.

`Fresa` is in the *Julia General Registry* <sup><a id="fnr.1" class="footref" href="#fn.1" role="doc-backlink">1</a></sup> and so installation is straightforward, with no dependencies beyond a small number of standard Julia packages.

\bibliographystyle{unsrtnat}
\bibliography{}


# Footnotes

<sup><a id="fn.1" href="#fnr.1">1</a></sup> <https://github.com/JuliaRegistries/General>
