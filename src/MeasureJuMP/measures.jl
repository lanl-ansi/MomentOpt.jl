"""
    AbstractGMPMeasure

Abstract type to represent measures in MeasureJuMP. Implementation of subtypes should implement MOI.get functions for `::PolynomialVariables`, `::Support`, `::RelaxationType`, and `::MomentBasis`.
"""
abstract type AbstractGMPMeasure end

MOI.get(m::AbstractGMPMeasure, ::PolynomialVariables) = m.variables
MOI.get(m::AbstractGMPMeasure, ::Support) = m.support
MOI.get(m::AbstractGMPMeasure, ::RelaxationType) = m.relax_type
MOI.get(m::AbstractGMPMeasure, ::MomentBasis) = m.moment_basis
MOI.get(m::AbstractGMPMeasure, ::MomentFunction) = nothing

function Base.:(==)(m1::AbstractGMPMeasure, m2::AbstractGMPMeasure)
    return all( get(m1, attr) == get(m2, attr) for attr in [PolynomialVariales(), Support(), RelaxationType(), MomentBasis()])
end

"""
    maxdegree_monomials(m::AbstractMeasure, max_degree::Int) 

Returns all monomials in the moment basis of m up to degree max_degree.    
"""
function maxdegree_monomials(m::AbstractMeasure, max_degree::Int) 
    return MB.maxdegree_basis(get.(m, [MomentBasis(), PolynomialVariables()])..., max_degree)  
end

"""
    covering_monomials(m::AbstractMeasure, p::MP.AbstractPolynomialLike) 

Returns all monomials in the moment basis of m covering the monomials of p.    
"""
function covering_monomials(m::AbstractMeasure,  p::MP.AbstractPolynomialLike) 
    return MB.basis_covering_monomials(get.(m, [MomentBasis(), PolynomialVariables()])..., monomials(p))  
end

"""
    moment_vector(basis::MB.AbstractPolynomialBasis, m::AbstractMeasure)

Returns the integral with respect to m over the polynomials in basis.
"""
function moment_vector(basis::MB.AbstractPolynomialBasis, m::AbstractMeasure)
    @assert get(m, ::RelaxationType()) isa EXACT_RELAXATION "Some requested moments are not available"
    return m.moment_function(basis)
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
    relax_type::NO_RELAXATION
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
    relax_type::EXACT_RELAXATION
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
    relax_type::EXACT_RELAXATION
    moment_basis::Nothing
    moment_function::Function
end

ZeroMeasure() = ZeroMeasure(nothing, nothing, EXACT_RELAXATION(), nothing, x -> 0)

export VariableMeasure

"""
    VariableMeasure{S, V, T}

Type representing a measure that can be relaxed via a conic relaxation. 
"""
struct VariableMeasure{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPMeasure
    support::S
    variables::Vector{V}
    relax_type::CONIC_FORMULATION
    moment_basis::T
end
