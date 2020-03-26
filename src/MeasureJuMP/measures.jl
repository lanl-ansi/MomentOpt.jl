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
abstract type AbstractGMPMeasure <: AbstractGMPVariable end

struct MeasureRef <: GMPVariableRef
    model::GMPModel
    index::MOI.VariableIndex
end

gmp_variable_refererence_type(m::AbstractGMPMeasure) = MeasureRef

MOI.get(m::AbstractGMPMeasure, ::Variables) = m.variables
MOI.get(m::AbstractGMPMeasure, ::Support) = m.support
MOI.get(m::AbstractGMPMeasure, ::ApproximationType) = m.relax_type
MOI.get(m::AbstractGMPMeasure, ::GMPBasis) = m.moment_basis

function MOI.get(m::AbstractGMPMeasure, ::MomentFunction) 
    get(m, ApproximationType()) == EXACT_APPROXIMATION ? m.moment_function : nothing
end

function Base.:(==)(m1::AbstractGMPMeasure, m2::AbstractGMPMeasure)
    return all( get(m1, attr) == get(m2, attr) for attr in [GenericMeasureAttributes..., MomentFunction()])
end

"""
    covering_basis(t::AbstractGMPVariable, p::MP.AbstractPolynomialLike)
    covering_basis(t::AbstractGMPVariable, p::Number)

Returns all monomials in the basis of t covering the monomials of p.    
"""
function covering_basis(t::AbstractGMPVariable, p::MP.AbstractPolynomialLike)
    return MB.basis_covering_monomials(get.(t, [GMPBasis(), Variables()])..., monomials(p)) 
end

"""
    mono_one(vars)

returns the monomial 1 defined on vars. 
"""
mono_one(vars:::Vector{MP.AbstractVariable}) = prod(var^0 for var in vars)
mono_one(::Nothing) = 1

covering_basis(t::AbstractGMPVariable, p::Number) = covering_basis(t, p*mono_one(get(t, Variables())))

"""
    maxdegree_basis(t::AbstractGMPVariable, d::Int)

Returns all monomials up to degree d in the basis of t.    
"""
function MB.maxdegree_basis(t::AbstractGMPVariable, d::Int)
    return maxdegree_basis(get.(t, [GMPBasis(), Variables()])..., d) 
end

"""
    moment_vector(basis::MB.AbstractPolynomialBasis, m::AbstractMeasure)

Returns the integral with respect to m over the polynomials in basis.
"""
function moment_vector(basis::MB.AbstractPolynomialBasis, m::AbstractMeasure)
    @assert get(m, ::ApproximationType()) isa EXACT_APPROXIMATION "Some requested moments are not available"
    return get(m, MomentFunction())(basis)
end

export integrate

"""
     integrate(p::Number, m::AbstractMeasure)
     integrate(p::MP.AbstractPolynomialLike, m::AbstractMeasure)

Returns the integral of p with respect to m. 
"""
function integrate(p::Number, m::AbstractMeasure)
    return p*moment_vector(maxdegree_monomials(m, 0))
end

function integrate(p::MP.AbstractPolynomialLike, m::AbstractMeasure)
    return dot(coefficients(p), moment_vector(covering_monomials(m, monomials(p))))
end

export SymbolicMeasure

"""
    SymbolicMeasure{S, V, T} <: AbstractGMPMeasure

Type representing a symbolic measure. This type does not allow to compute integrals or relaxations.
"""
struct SymbolicMeasure{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPMeasure
    support::S
    variables::Vector{V}
    relax_type::NO_APPROXIMATION
    moment_basis::T
end

function SymbolicMeasure(support, variables, moment_basis) 
    return SymbolicMeasure(support, variables, NO_RELAXATION(), moment_basis) 
end

export AnalyticMeasure

"""
    AnalyticMeasure{V, T} <: AbstractGMPMeasure

Type representing an analytiv measure. Its field moment_function allows to integrate polynomials. 
"""
struct AnalyticMeasure{V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} <: AbstractGMPMeasure
    support::Nothing
    variables::Vector{V}
    relax_type::EXACT_APPROXIMATION
    moment_basis::T
    moment_function::Function
end

function AnalyticMeasure(variables, moment_basis, moment_function) 
    return AnalyticMeasure(nothing, variables, EXACT_RELAXATION(), moment_basis, moment_function) 
end

export ZeroMeasure

"""
    ZeroMeasure <: AbstractGMPMeasure
   
Type representing the measure giving zero weight to anything. 
"""
struct ZeroMeasure <: AbstractGMPMeasure
    support::Nothing
    variables::Nothing
    relax_type::EXACT_APPROXIMATION
    moment_basis::Nothing
    moment_function::Function
end

ZeroMeasure() = ZeroMeasure(nothing, nothing, EXACT_RELAXATION(), nothing, x -> 0)

function covering_basis(t::ZeroMeasure, p::MP.AbstractPolynomialLike)
    return MB.basis_covering_monomials(get.(t, MonomialBasis, variables(p)), monomials(p)) 
end

export VariableMeasure

"""
    VariableMeasure{S, V, R, T}

Type representing a measure that can be relaxed via a conic relaxation. 
"""
struct VariableMeasure{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, R <: AbstractApproximation, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPMeasure
    support::S
    variables::Vector{V}
    relax_type::R
    moment_basis::T
end



