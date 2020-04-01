Base.broadcastable(o::AbstractGMPObject) = Ref(o)

"""
    AbstractApproximationType

Defines what kind of certificate or relaxation is used for the approximation of a measure or a function.
"""
abstract type AbstractApproximationType end
struct NO_APPROXIMATION <: AbstractApproximationType end
struct EXACT_APPROXIMATION <: AbstractApproximationType end
abstract type AbstractApproximation <: AbstractApproximationType end
struct DefaultApproximation <: AbstractApproximation end #TODO to be removed as soon as default approximation is available.

# GMPObjectAttributes and get functions

"""
    Variables()

Attribute for the variables a measure is acting on. 
"""
struct Variables <: MOI.AbstractVariableAttribute end

"""
    BSASet()

Attribute for the basic semialgebraic set related to a GMPObject.
"""
struct BSASet <: MOI.AbstractVariableAttribute end

"""
    ApproximationType()

Attribute for the type of approximation used for a measure/continuous function.
"""
struct ApproximationType <: MOI.AbstractVariableAttribute end

"""
    ApproximationBasis()

Attribute for the type of MultivariateBases, used to express the moment of a measure or the monomials of a polynomial.
"""
struct ApproximationBasis <: MOI.AbstractVariableAttribute end

"""
    ApproximationFunction()

Attribute for a function allowing to approximate the action of a GMPObject.
"""
struct ApproximationFunction <: MOI.AbstractVariableAttribute end

const GenericObjectAttributes = [Variables(), BSASet(), ApproximationType(), ApproximationBasis()]

MOI.get(m::AbstractGMPObject, ::Variables) = m.variables
MOI.get(m::AbstractGMPObject, ::BSASet) = m.bsa_set
MOI.get(m::AbstractGMPObject, ::ApproximationType) = m.approx_type
MOI.get(m::AbstractGMPObject, ::ApproximationBasis) = m.approx_basis
function MOI.get(m::AbstractGMPObject, ::ApproximationFunction) 
    get(m, ApproximationType()) == EXACT_APPROXIMATION ? m.approx_func : nothing
end

function Base.:(==)(m1::AbstractGMPObject, m2::AbstractGMPObject)
    return typeof(m1) == typeof(m2) && all(MOI.get(m1, attr) == MOI.get(m2, attr) for attr in GenericObjectAttributes)
end

# GMPObjectAttributes functionalities

"""
    covering_basis(t::AbstractGMPObject, p::MP.AbstractPolynomialLike)
    covering_basis(t::AbstractGMPObject, p::Number)

Returns all monomials in the basis of t covering the monomials of p.    
"""
function covering_basis(t::AbstractGMPObject, p::MP.AbstractPolynomialLike)
    return MB.basis_covering_monomials(get.(t, [GMPBasis(), Variables()])..., monomials(p)) 
end

"""
    _mono_one(vars)

returns the monomial 1 defined on vars. 
"""
_mono_one(vars::Vector{MP.AbstractVariable}) = prod(var^0 for var in vars)
_mono_one(::Nothing) = 1

covering_basis(t::AbstractGMPObject, p::Number) = covering_basis(t, p*mono_one(get(t, Variables())))


"""
    maxdegree_basis(t::AbstractGMPObject, d::Int)

Returns all monomials up to degree d in the basis of t.    
"""
function MB.maxdegree_basis(t::AbstractGMPObject, d::Int)
    return maxdegree_basis(get.(t, [GMPBasis(), Variables()])..., d) 
end

"""
    eval_vector(basis::MB.AbstractPolynomialBasis, m::AbstractObject)

Returns the values of m corresponding to the monomials in basis. 
"""
function eval_vector(basis::MB.AbstractPolynomialBasis, m::AbstractGMPObject)
    @assert get(m, ApproximationType()) isa EXACT_APPROXIMATION
    return get(m, ApproximationFunction()).(basis)
end

"""
    AbstractGMPMeasure

Abstract type to represent measures in MeasureJuMP. Implementation of subtypes should implement MOI.get functions for `::Variables`, `::Support`, `::ApproximationType`, and `::GMPBasis`.
if get(m, ApproximationType()) == EXACT_APPROXIMATION, the measure m should have a field m.moment_function::Function. 
By default there are four implementations of AbstractGMPMeasure:
  * [SymbolicMeasure](@ref) :
  * [AnalyticMeasure](@ref) :
  * [ZeroMeasure](@ref) :
  * [VariableMeasure](@ref) :

"""
abstract type AbstractGMPMeasure <: AbstractGMPObject end
include("MeasureJuMP/measures.jl")

"""
    AbstractGMPContinuous

Under some regularity assumptions one can say that the dual space of measures is the space of continuous functions. This is for example the case if the support of the measure is compact. This abstract type allows to define custom testfunctions.
"""
abstract type AbstractGMPContinuous <: AbstractGMPObject end
include("MeasureJuMP/continuous.jl")
