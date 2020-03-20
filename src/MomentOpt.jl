module MomentOpt

import Reexport

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
Reexport.@reexport using JuMP
const PJ = PolyJuMP

Reexport.@reexport using SumOfSquares

using OrderedCollections
using LinearAlgebra

const MT = Union{Number, AbstractPolynomialLike}



include("measurevariable.jl")
include("knownmeasures.jl")

abstract type AbstractMeasureRef <: JuMP.AbstractVariableRef end
include("measexpr.jl")
include("momexpr.jl")

abstract type AbstractGMPConstraint <: JuMP.AbstractConstraint end
include("momentconstraints.jl")
include("measureconstraints.jl")
include("substitutionconstraints.jl")

abstract type AbstractGMPConstraintRef end
include("gmpmodel.jl")

include("measureref.jl")
include("constraintref.jl")

include("measuremacro.jl")
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
