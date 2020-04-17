module MomentOpt

using Reexport

using MultivariatePolynomials
const MP = MultivariatePolynomials
Reexport.@reexport using MultivariateBases
const MB = MultivariateBases
include("MBext/MBextra.jl")
include("MPext/MPextra.jl")

# MOI extension
using MathOptInterface
const MOI = MathOptInterface

using MutableArithmetics
const MA = MutableArithmetics
Reexport.@reexport using JuMP

using LinearAlgebra

Reexport.@reexport using SemialgebraicSets
include("SASext/SASextra.jl")

Reexport.@reexport using SumOfSquares


# begin MomentOpt
include("approximation.jl")
include("approximationscheme.jl")
include("objects.jl")
include("defaultmeasures.jl")
include("variables.jl")
include("gmpaffexpr.jl")
include("affexpr.jl")

abstract type GMPSubstitution <: JuMP.AbstractJuMPScalar end
# define MomentSubstitution
# define IntegralSubstitution

include("constraints.jl")

#
include("gmpmodel.jl")

include("approximate.jl")

include("MMext/MMextra.jl")

include("gmppostproc.jl")

end# module
