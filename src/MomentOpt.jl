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

abstract type AbstractGMPSet <: MOI.AbstractScalarSet end
abstract type AbstractGMPShape <: JuMP.AbstractShape end
abstract type AbstractGMPConstraint <: JuMP.AbstractConstraint end


abstract type GMPSubstitution <: JuMP.AbstractJuMPScalar end
# define MomentSubstitution
# define IntegralSubstitution

include("constraints.jl")
#
include("gmpmodel.jl")

Reexport.@reexport using SumOfSquares

include("approximate.jl")

include("MMext/MMextra.jl")
#include("relax.jl")

import LinearAlgebra.eigen
include("gmppostproc.jl")


end# module
