"""
    AbstractGMPContinuous

Under some regularity assumptions one can say that the dual space of measures is the space of continuous functions. This is for example the case if the support of the measure is compact. This abstract type allows to define custom testfunctions.
"""
abstract type AbstractGMPContinuous <: AbstractGMPVariable end

struct ContinuousRef <: GMPVariableRef
    model::GMPModel
    index::MOI.VariableIndex
end

gmp_variable_refererence_type(m::AbstractGMPContinuous) = ContinuousRef


MOI.get(m::AbstractGMPContinuous, ::Variables) = m.variables
MOI.get(m::AbstractGMPContinuous, ::Domain) = m.domain
MOI.get(m::AbstractGMPContinuous, ::ApproximationType) = m.stregthen_type
MOI.get(m::AbstractGMPContinuous, ::GMPBasis) = m.monom_basis
function MOI.get(m::AbstractGMPMeasure, ::CoefFunction) 
    get(m, ApproximationType()) == EXACT_APPROXIMATION ? m.coef_function : nothing
end

function Base.:(==)(m1::AbstractGMPContinuous, m2::AbstractGMPContinuous)
    return all( get(m1, attr) == get(m2, attr) for attr in [GenericContinuousAttributes..., CoefFunction()])
end

"""
    coefficient_vector(basis::MB.AbstractPolynomialBasis, f::AbstractContinuous)

Returns the coefficients with respect to the polynomials in basis, approximating f.
"""
function coefficient_vector(basis::MB.AbstractPolynomialBasis, f::AbstractContinuous)
    @assert get(f, ::ApproximationType()) isa EXACT_APPROXIMATION "Some requested coefficients are not available"
    return get(f, CoefFunction())(basis)
end

export approximate
"""
    approximate(f::AbstractContinuous, max_degree::Int)

Returns a polynomial approximation of f of degree max_degree.
"""
function approximate(f::AbstractContinuous, max_degree::Int)
    basis = MB.maxdegree_basis(f, max_degree)
    return dot(coefficient_vector(basis), basis)
end


export SymbolicContinuous

"""
    SymbolicContinuous{S, V, T} <: AbstractGMPContinuous

Type representing a symbolic continuous function. This type does not allow to compute integrals or relaxations.
"""
struct SymbolicContinuous{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    domain::S
    variables::Vector{V}
    stregthen_type::NO_APPROXIMATION
    monom_basis::T
end

function SymbolicContinuous(support, variables, monom_basis) 
    return SymbolicContinuous(support, variables, NO_RELAXATION(), monom_basis) 
end

export AnalyticContinuous

"""
    AnalyticContinuous{V, T} <: AbstractGMPContinuous

Type representing an analytic continuous function. Its field coef_function allows to compute the coefficients for arbitrary close polynomial approximations.
"""
struct AnalyticContinuous{V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    domain::Nothing
    variables::Vector{V}
    stregthen_type::EXACT_RELAXATION
    monom_basis::T
    coef_function::Function
end

function AnalyticContinuous(variables, monom_basis, monom_function) 
    return AnalyticContinuous(nothing, variables, EXACT_RELAXATION(), monom_basis, monom_function) 
end

export ConstantContinuous

"""
    ConstantContinuous <: AbstractGMPContinuous
   
Type representing constant functions  
"""
struct ConstantContinuous <: AbstractGMPContinuous
    domain::Nothing
    variables::Nothing
    stregthen_type::EXACT_RELAXATION
    monom_basis::Nothing
    coef_function::Function
end

ConstantContinuous(a::Number) = ConstantContinuous(nothing, nothing, EXACT_RELAXATION(), nothing, x -> x == 1 ? a : 0)
ZeroContinuous() = ConstantContinuous(0)
OneContinuous() = ConstantContinuous(1)

function approximate(f::ConstantContinuous, ::Int)
    return get(f, CoefFunction())(1)*mono_one(get(f, Variables()))
end

export VariableContinuous

"""
    VariableContinuous{S, V, R, T}

Type representing a continuous funciton that can be stregthend via a conic approximation. 
"""
struct VariableContinuous{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, R <: AbstractApproximation, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    domain::S
    variables::Vector{V}
    stregthen_type::R
    monom_basis::T
end


