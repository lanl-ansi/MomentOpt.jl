module MomentOpt

using Reexport

using MathOptInterface
const MOI = MathOptInterface

using MutableArithmetics
const MA = MutableArithmetics

using LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials
include("MPext/MPextra.jl")

Reexport.@reexport using JuMP

Reexport.@reexport using SemialgebraicSets
include("SASext/SASextra.jl")

using SumOfSquares

Reexport.@reexport using MultivariateMoments
const MM = MultivariateMoments

Reexport.@reexport using MultivariateBases
const MB = MultivariateBases
include("MBext/MBextra.jl")

abstract type AbstractGMPModel <: JuMP.AbstractModel end


include("approximation.jl")
include("objects.jl")
include("variables.jl")
include("defaultmeasures.jl")
include("gmpaffexpr.jl")
include("affexpr.jl")
include("approximationscheme.jl")

abstract type GMPSubstitution <: JuMP.AbstractJuMPScalar end
# define MomentSubstitution

include("constraints.jl")
#
include("gmpmodel.jl")
include("MMext/MMextra.jl")
include("approximate.jl")
include("gmppostproc.jl")

end# module
