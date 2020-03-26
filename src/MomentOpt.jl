module MomentOpt

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases

# MOI extension

using MathOptInterface
const MOI = MathOptInterface
include("MeasureJuMP/attributes.jl")

using Reexport

Reexport.@reexport using JuMP

abstract type AbstractGMPVariable  <: JuMP.AbstractVariable end
const GMPAffExpr{T <: Number, V <: AbstractGMPVariable} = JuMP.GenericAffExpr{T, V}
abstract type AbstractGMPConstraint <: AbstractConstraint end

include("gmpmodel.jl")

#=

struct NoScalar <: AbstractGMPScalar end
function JuMP.function_string(io::IO, s::AbstractGMPScalar) end
function degree(s::AbstractGMPScalar) end 

abstract type AbstractGMPScalarSet <: JuMP.AbstractScalarSet 
include("MeasureJuMP/sets.jl")


abstract type AbstractGMPConstraint end <: AbstractConstraint
function constraint_function(c::AbstractGMPConstraint) end
function JuMP.shape(x::AbstractGMPConstraint) end
function JuMP.constraint_string(print_mode, constraint::AbstractGMPConstraint) end


import LinearAlgebra.dot
using SemialgebraicSets
include("variables.jl")



Reexport.@reexport using SumOfSquares
=#
#=
using MutableArithmetics
const MA = MutableArithmetics

using MultivariatePolynomials
const MP = MultivariatePolynomials

using MultivariateBases
const MB = MultivariateBases

using MultivariateMoments
const MM = MultivariateMoments

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

using PolyJuMP
const PJ = PolyJuMP

Reexport.@reexport using SumOfSquares

using OrderedCollections
using LinearAlgebra

const MT = Union{Number, MP.AbstractPolynomialLike}


include("info.jl")
include("gmpmodel.jl")

include("measurevariable.jl")
include("knownmeasures.jl")

include("measureexpr.jl")
include("momentexpr.jl")

include("momentobjective.jl")

abstract type AbstractGMPConstraint <: JuMP.AbstractConstraint end
include("momentconstraints.jl")
include("measureconstraints.jl")
include("substitutionconstraints.jl")

include("measureref.jl")
include("constraintref.jl")

include("measuremacro.jl")
=#
#=
include("momcon.jl")

include("momentsequence.jl")
include("riesz.jl")
include("relaxationtypes.jl")
include("meas.jl")
include("model.jl")
include("relax.jl")
include("postproc.jl")
include("macros.jl")
=#
end# module
