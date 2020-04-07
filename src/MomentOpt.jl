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
include("defaultmeasures.jl")

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

include("constraints.jl")
#
include("gmpmodel.jl")


Reexport.@reexport using SumOfSquares

include("MPextra.jl")
include("approximate.jl")
#include("relax.jl")

#=
include("postproc.jl")
=#

end# module
