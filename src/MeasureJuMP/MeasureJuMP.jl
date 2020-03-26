module MeasureJuMP

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases

include("formulationtypes.jl")

# MOI extension

using MathOptInterface
const MOI = MathOptInterface
include("attributes.jl")

import LinearAlgebra.dot
using SemialgebraicSets

"""
    mono_one(vars)

returns the monomial 1 defined on vars. 
"""
mono_one(vars:::Vector{MP.AbstractVariable}) = prod(var^0 for var in vars)
mono_one(::Nothing) = 1

abstract type AbstractGMPType end
include("measures.jl")
include("continuous.jl")

include("sets.jl")
# JumP extension

using JuMP

include("variable.jl")

end # module
