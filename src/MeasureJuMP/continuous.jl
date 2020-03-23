"""
    AbstractGMPContinuous

Under some regularity assumptions one can say that the dual space of measures is the space of continuous functions. This is for example the case if the support of the measure is compact. This abstract type allows to define custom testfunctions.
"""
abstract type AbstractGMPContinuous end

#TODO introduce correct names for attributes
#
MOI.get(m::AbstractGMPContinuous, ::PolynomialVariables) = m.variables
MOI.get(m::AbstractGMPContinuous, ::Domain) = m.domain
MOI.get(m::AbstractGMPContinuous, ::RelaxationType) = m.relax_type
MOI.get(m::AbstractGMPContinuous, ::MomentBasis) = m.moment_basis
MOI.get(m::AbstractGMPContinuous, ::MomentFunction) = nothing

function Base.:(==)(m1::AbstractGMPContinuous, m2::AbstractGMPContinuous)
    return all( get(m1, attr) == get(m2, attr) for attr in [PolynomialVariales(), Domain(), RelaxationType(), MomentBasis(), MomentFuntion()])
end

export SymbolicContinuous

"""
    SymbolicContinuous{S, V, T} <: AbstractGMPContinuous

Type representing a symbolic continuous function. This type does not allow to compute integrals or relaxations.
"""
struct SymbolicContinuous{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    domain::S
    variables::Vector{V}
    relax_type::NO_RELAXATION
    moment_basis::T
end

function SymbolicContinuous(support, variables, moment_basis) 
    return SymbolicContinuous(support, variables, NO_RELAXATION(), moment_basis) 
end

export AnalyticContinuous

"""
    AnalyticContinuous{V, T} <: AbstractGMPContinuous

Type representing an analytic continuous function. Its field moment_function allows to evluate polynomials. 
"""
struct AnalyticContinuous{V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    domain::Nothing
    variables::Vector{V}
    relax_type::EXACT_RELAXATION
    moment_basis::T
    moment_function::Function
end

function AnalyticContinuous(variables, moment_basis, moment_function) 
    return AnalyticContinuous(nothing, variables, EXACT_RELAXATION(), moment_basis, moment_function) 
end

export ConstantContinuous

"""
    ConstantContinuous <: AbstractGMPContinuous
   
Type representing constant functions  
"""
struct ConstantContinuous <: AbstractGMPContinuous
    domain::Nothing
    variables::Nothing
    relax_type::EXACT_RELAXATION
    moment_basis::Nothing
    moment_function::Function
end

ConstantContinuous(a::Number) = ConstantContinuous(nothing, nothing, EXACT_RELAXATION(), nothing, x -> a)
ZeroContinuous() = ConstantContinuous(0)
OneContinuous() = ConstantContinuous(1)

export VariableContinuous

"""
    VariableContinuous{S, V, T}

Type representing a measure that can be relaxed via a conic relaxation. 
"""
struct VariableContinuous{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    domain::S
    variables::Vector{V}
    relax_type::CONIC_FORMULATION
    moment_basis::T
end


