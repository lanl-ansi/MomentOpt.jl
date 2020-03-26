# MomentOpt.jl - Modelization and Conic Relaxations for Generalized Moment Problems

MomentOpt.jl is a Julia package to model Generalized Moment Problems and to approximate solutions via conic relaxations. 

The two main ideas of MomentOpt.jl are 
  * to enable users to easily formulate their problem in terms of measures and
  * to compute approximations based on different conic relaxations.

MomentOpt.jl uses the MultivariatePolynomials.jl interface to represent polynomials and moments. We recommend using the implementation DynamicPolynomials.jl. Polynomial variables are then created with the `DynamicPolynomials.@polyvar` macro

```julia
using MomentOpt
using DynamicPolynomials

@polyvar x y
```

MomentOpt.jl is implemented as a JuMP.jl extension. In particular it uses the same syntax for modeliztion. For example:
  * `m = GMPModel()` generates an empty model representing a generalized moment problem
  * `@measure m [x, y] args...` adds a measure variable to m. This macro behaves like the `JuMP.@variable` macro.  
  * `@constraint m args...` adds constraints to m. This macro is actually the macro `JuMP.@constraint`.

```@contents
    Pages = ["gmp.md"]
    Depth = 2
```
