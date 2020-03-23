module MomentOpt

using Reexport

Reexport.@reexport using JuMP
include("info.jl")
include("gmpmodel.jl")


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
