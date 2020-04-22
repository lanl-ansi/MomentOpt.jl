module MomentOpt

using Reexport
# MOI extension
using MathOptInterface
const MOI = MathOptInterface

using MutableArithmetics
const MA = MutableArithmetics
Reexport.@reexport using JuMP

using LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials
Reexport.@reexport using SemialgebraicSets
include("SASext/SASextra.jl")

Reexport.@reexport using SumOfSquares
# extending MultivariateMoments.jl for the use of MomentOpt.jl
Reexport.@reexport using MultivariateMoments
const MM = MultivariateMoments



Reexport.@reexport using MultivariateBases
const MB = MultivariateBases
include("MBext/MBextra.jl")
include("MPext/MPextra.jl")

include("approximation.jl")

include("objects.jl")
include("defaultmeasures.jl")
include("variables.jl")
include("gmpaffexpr.jl")
include("affexpr.jl")



# begin MomentOpt
include("approximationscheme.jl")

abstract type GMPSubstitution <: JuMP.AbstractJuMPScalar end
# define MomentSubstitution
# define IntegralSubstitution

include("constraints.jl")

#
include("gmpmodel.jl")

include("MMext/MMextra.jl")

include("approximate.jl")



include("gmppostproc.jl")

end# module
