module MomentOpt

using Reexport

using MultivariatePolynomials
const MP = MultivariatePolynomials
Reexport.@reexport using MultivariateBases
const MB = MultivariateBases

# MOI extension
using MathOptInterface
const MOI = MathOptInterface
include("MeasureJuMP/attributes.jl") # maybe move to gmpmodel

Reexport.@reexport using JuMP

import LinearAlgebra.dot
Reexport.@reexport using SemialgebraicSets

abstract type AbstractGMPObject end
include("objects.jl")

abstract type AbstractGMPVariable  <: JuMP.AbstractVariable end
abstract type AbstractGMPVariableRef <: JuMP.AbstractVariableRef end
include("variables.jl")

abstract type AbstractGMPScalar <: JuMP.AbstractJuMPScalar end
include("gmpaffexpr.jl")
include("affexpr.jl")

abstract type GMPSubstitution <: JuMP.AbstractJuMPScalar end
# define MomentSubstitution
# define IntegralSubstitution

abstract type AbstractGMPSet <: MOI.AbstractScalarSet end
#include("MeasureJuMP/sets.jl")

abstract type AbstractGMPShape <: JuMP.AbstractShape end
#include("shapes.jl")

abstract type AbstractGMPConstraint <: JuMP.AbstractConstraint end
# define LinearMeasureConstraint
# define LinearContinuousConstraint
# define LinearMomentConstraint
# define LinearIntegralConstriant
# define MomentSubstitutionConstraint
# define IntegralSubstitutionConstraint
include("constraints.jl")
#
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
