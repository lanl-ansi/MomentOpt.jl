# Getting started

## Installation

MomentOpt is a package for [Julia](https://julialang.org/). In order to use it you need to download and install Julia first. Instructions on how to do this can be found [here](https://julialang.org/downloads/). 

Once Julia is installed, MomentOpt can be added through Julia's internal package manager.
```julia
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.4.1 (2020-04-1)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> ]

pkg> add MomentOpt
```

MomentOpt uses [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) to model polynomial variables. The user can choose between different implementations for this polynomial interface. We recommend DynamicPolynomials:

```julia
pkg> add DynamicPolynomials
```

Having installed these two packages we could start modelling with MomentOpt. However, to be able to compute approximations, we need to install at least one solver. A comprehensive list of possible solvers can be found on [the JuMP website](https://github.com/JuliaOpt/JuMP.jl). 
As an open source SDP solver we could install CSDP:
```julia
pkg> add CSDP
```

Once these packages have been added, they do not need to be added again when starting a new session. In case you are using [Jupyter notebooks](https://jupyter.org/) to run Julia, the syntax might slightly vary. Please refer to the [IJulia](https://github.com/JuliaLang/IJulia.jl) documentation. 


## First steps
To start your session open julia and load MomentOpt, an implementation of MultivariatePolynomials and an SDP solver, e.g.
```julia
julia> using MomentOpt
julia> using DynamicPolynomials
julia> using CSDP
```

The first thing you would probably want to do, is to define the variables your measures will be acting on. In DynamicPolynomials we can use the `@polyvar` macro to do so:
```julia
julia> @polyvar x y z
(x, y, z)
```

Now we can start modeling our generalized moment problem. The modeling process is consitent with [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl). In particular, we always start by building a model:
```julia
julia> gmp = GMPModel()
A JuMP Model
Feasibility problem with:
Variables: 0
Constraints: 0
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 0
Degree for approximation 0
Solver for approximation: No optimizer attached.
```

As you can see gmp is an abstract JuMP model that currently represents a feasibility problem without any variables or constraints. The next three lines consider the specifications for the approximation and will be discussed later.
Note that no optimizer has been attached to the model. To do so, we can do
```julia
julia> set_optimizer(gmp, CSDP.Optimizer)
julia> gmp
A JuMP Model
Feasibility problem with:
Variables: 0
Constraints: 0
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 0
Degree for approximation 0
Solver for approximation: CSDP
```
Alternatively, we could have attached the optimizer from the beginning by defining
```julia
julia> gmp = GMPModel(CSDP.Optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Constraints: 0
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 0
Degree for approximation 0
Solver for approximation: CSDP
```
The variables of a generalized moment problem are measures. To add a measure to our model, we can do
```julia
julia> @variable gmp mu Meas([x, y])
mu
```
This is following the JuMP syntax, i.e., we first specify the model, than a name for our variable, and then specify that the new variable represents a measure acting on the polynomial variables x and y. It is mandatory to define on which variables a measure is acting. The polynomial variables of a measure variable can be querried via 
```julia
julia> variables(mu)
2-element Array{PolyVar{true},1}:
 x
 y
```
The most important keyword argument of the `Meas()` constructor is `support` which determines the semialgebraic set which contains the support of the measure variable:
```julia
julia> @variable gmp nu Meas([x,y], support = @set(1-x^2-y^2>=0))
nu
julia> support(nu)
Basic semialgebraic Set defined by no equality
1 inequalitty
 -x^2 - y^2 + 1 ≥ 0
```
For more information on generating semi-algebraic sets via the `@set` macro we refer to [SemiAlgebraicSets.jl](https://github.com/JuliaAlgebra/SemialgebraicSets.jl).
The remaining keywords are discussed in the documentation of [`Meas`](@ref). 

Note that you can create anonymous measure variables or containers of measure variables
```julia
julia> phi = @variable gmp [1:2] Meas([x,y,z])
2-element Array{GMPVariableRef{MomentOpt.AbstractGMPMeasure},1}:
 noname
 noname
julia> vars = [x, y, z]
julia> psi = @variable gmp [i = 1:3] Meas([vars[i]]; support = @set(vars[i]^2 <= 1))
3-element Array{GMPVariableRef{MomentOpt.AbstractGMPMeasure},1}:
 noname
 noname
 noname
julia> support.(psi)
3-element Array{BasicSemialgebraicSet{Int64,Polynomial{true,Int64},FullSpace},1}:
 { (x) | -x^2 + 1 ≥ 0 }
 { (y) | -y^2 + 1 ≥ 0 }
 { (z) | -z^2 + 1 ≥ 0 }
```
Note that the `.` behind `support` broadcasts the function to each element of the vector `psi`. 
If we look at our model again, we find that mu and nu have been registered in the model.
```julia 
julia> gmp
A JuMP Model
Feasibility problem with:
Variables: 7
Constraints: 0
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 2
Degree for approximation 2
Solver for approximation: CSDP
Names registered in the model: mu, nu
```
We cannot add a new object with names that are alrady registered:
```
julia> @variable gmp mu Meas([x, y])
ERROR: An object of name mu is already attached to this model. If this is intended, 
consider using the anonymous construction syntax, e.g., x = @variable(model, [1:N], ...) 
where the name of the object does not appear inside the macro.
```

 Also note, that we have registered 6 variables and that the maximum degree of data is 2 now, since the definition of the measure `nu` involved a polynomial of degree 2. The maximum degree of data will update automatically, when adding variables or constraints. The degree for approximation is always at least as large as the degree of the data and is updated automatically:
```julia
julia> set_approximation_degree(gmp, 4)
julia> gmp
A JuMP Model
Feasibility problem with:
Variables: 7
Constraints: 0
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 2
Degree for approximation 4
Solver for approximation: CSDP
Names registered in the model: mu, nu

julia> set_approximation_degree(gmp, 1)
┌ Warning: Requested approximation degree 1 is too low to cover all data. 
The approximation degree has been set to the minimal value possible and now is 2.
```

Next we can add constraints to our model. However, we will start with a new model. A typical constraint for generalized moment problems is that the measure variable is a probability measure, which is added as follwos
```julia
julia> gmp = GMPModel(CSDP.Optimizer);
julia> @variable gmp mu Meas([x], support = @set(x-x^2>=0));
julia> @constraint gmp Mom(1, mu) == 1
⟨1, mu⟩ = 1.0
```
Here `Mom(1, mu)` expressed the integral of `1` with respect to `mu`. Constraints can have names
```julia
julia> @constraint gmp xcon Mom(x, mu) == 0.5;
julia> xcon
xcon : ⟨x, mu⟩ = 0.5
julia> @constraint gmp xicon[i = 2:4] Mom(x^i, mu) .== 1/(i+1);
julia> xicon[2]
⟨x², mu⟩ = 0.3333333333333333

```

So far we have only considered moment constraints. MomentOpt also supports measure constraints:
```julia
julia> lambda = lebesgue_measure_box(x => [0, 1])
AnalyticMeasure
```
defines an AnalyticMeasure, representing the Lebesgue measure on [0, 1]. It can be used to integrate polynomials
```julia
julia> integrate(x^2 - x + 1, lambda)
0.8333333333333333
julia integrate.(monoimals([x], 0:4), lambda)
5-element Array{Float64,1}:
 0.2
 0.25
 0.3333333333333333
 0.5
 1.0
```
but it can also be used to define measure constraints:
```julia 
julia> gmp = GMPModel(CSDP.Optimizer);
julia> @variable gmp mu Meas([x,y], support = @set(1-x^2-y^2>=0));
julia> @variable gmp nu Meas([x,y], support = @set(1-x^2>= 0 && 1-y^2>=0));
julia> gmp
A JuMP Model
Feasibility problem with:
Variables: 2
Constraint: 1
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 2
Degree for approximation 2
Solver for approximation: CSDP
Names registered in the model: mu, nu
julia> lambda = lebesgue_measure_box(x => [-1, 1], y => [-1, 1]; normalize = true)
julia> mc = @constraint gmp mu + nu == lambda
mu + nu = AnalyticMeasure
julia> gmp
A JuMP Model
Feasibility problem with:
Variables: 2
Constraint: 1
Approxmation mode: DUAL_STRENGTHEN_MODE()
Maximum degree of data: 2
Degree for approximation 2
Solver for approximation: CSDP
Names registered in the model: mu, nu
```
Note that `mc` did not change the degree of data. The degree of a measure constraint is always equal to the approximation degree and determined at the optimization call.

To finish this introduction on modeling with MomentOpt, we have to set an objective function. This can be done via
```julia
julia> @objective gmp Max Mom(1, mu) 
⟨1, mu⟩
```

The model now can be approximated by setting an approximation order (it defaults to the maximum degree of data). In this example we also suppress the output of CSDP:
```julia
julia> set_approximation_degree(gmp, 10)
julia> set_silent(gmp);
julia> optimize!(gmp)
```
We can querry the stats of the solved approximation
```julia
julia> objective_value(gmp)
0.9369024582093924
julia> termination_status(gmp)
OPTIMAL::TerminationStatusCode = 1
julia> raw_status(gmp)
"Problem solved to optimality."
```
The solution can be expected by calling 
```julia
julia> measure(mu)
MultivariateMoments.Measure{Float64,Monomial{true},MonomialVector{true}}([0.062044293013069955, 1.818408257153911e-19, 0.01273398937697667, 4.4781213912861547e-20, 0.0064424590352839825, 2.0217159884897534e-21, 0.006442459035519992, -4.815750725967947e-20, 0.012733989378712668, -2.0895884615236935e-19  …  -3.8109683913485666e-18, 1.865988294254163e-18, -3.515461570762854e-18, 2.056500523116538e-18, 0.2816658397255027, -7.275476044869635e-21, 0.28166583972970993, -5.312847305918446e-18, 2.8599527338155522e-18, 0.9369024582093924], Monomial{true}[x¹⁰, x⁹y, x⁸y², x⁷y³, x⁶y⁴, x⁵y⁵, x⁴y⁶, x³y⁷, x²y⁸, xy⁹  …  x³, x²y, xy², y³, x², xy, y², x, y, 1])
```
This returns a sequence of moments. The dual to our constraint `mc` is a polynomial of degree 10 in the variables x and y.
```julia
julia> poly = dual(mc);
julia> maxdegree(poly)
10
juila> variables(poly) 
2-element Array{PolyVar{true},1}:
 x
 y
```

