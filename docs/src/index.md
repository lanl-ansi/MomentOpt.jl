# MomentOpt.jl - Modelization and Conic Relaxations for Generalized Moment Problems

MomentOpt.jl is a Julia package to model Generalized Moment Problems and to approximate solutions via conic relaxations as it is described in [Moments, Positive Polynomials and Their Applications](https://homepages.laas.fr/lasserre/book-flyer.pdf) or more recently in [The Moment-SOS Hierarchy](https://www.worldscientific.com/worldscibooks/10.1142/q0252).

The two main ideas of MomentOpt.jl are
  * use a high level syntax to define Generalized Moment Problems easily
  * provide different options to approximate solutions and switch between different formulations easily.

MomentOpt.jl is implemented as a [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) extension. In particular it uses the same syntax for modelization. For example:
 * `m = GMPModel()` generates an empty model representing a generalized moment problem
  * `@variable m Meas([x, y]; kwargs...) args...` adds a measure variable to `m`.
  * `@constraint m args...` adds constraints to `m`.
  * `set_optimizer(m, optimizer)` can be used to set the optimizer. Alternatively, one can define the optimizer from the beginning `m = GMPModel(optimizer)`.
  *  `optimize!(m)` is used to approximate a solution to `m`

As a JuMP extension all numerical solvers available through JuMP are available to be used in MomentOpt, too.

MomentOpt.jl uses the [MultivariatePolynomials.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) interface to represent polynomials and moments. We recommend using the implementation [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl).


# Content
```@contents
    Pages = ["gmp.md"]
    Depth = 2
```

# How to cite
See [citation.bib](https://www.github.com/lanl-ansi/MomentOpt.jl/blob/master/citation.bib).
