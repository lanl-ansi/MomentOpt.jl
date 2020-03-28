export approximate

"""
    approximate(f::AbstractContinuous, max_degree::Int)

Returns a polynomial approximation of f of degree max_degree.
"""
function approximate(f::AbstractContinuous, max_degree::Int)
    basis = MB.maxdegree_basis(f, max_degree)
    return dot(eval_vector(basis), basis)
end

export SymbolicContinuous

"""
    SymbolicContinuous{V, S, T} <: AbstractGMPContinuous

Type representing a symbolic continuous function. This type does not allow to compute integrals or relaxations.
"""
struct SymbolicContinuous{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::S
    approx_type::NO_APPROXIMATION
    approx_basis::T
end

function SymbolicContinuous(variables::Vector{V}, domain::S, monom_basis::T) where {S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} 
    return SymbolicContinuous(support, variables, NO_RELAXATION(), monom_basis) 
end

export AnalyticContinuous

"""
    AnalyticContinuous{V, T} <: AbstractGMPContinuous

Type representing an analytic continuous function. Its field coef_function allows to compute the coefficients for arbitrary close polynomial approximations.
"""
struct AnalyticContinuous{V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::Nothing
    approx_type::EXACT_RELAXATION
    approx_basis::T
    approx_function::Function
end

function AnalyticContinuous(variables::Vector{T}, monom_basis::T, monom_function::Function) where {V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} 
    return AnalyticContinuous(variables, nothing,  EXACT_RELAXATION(), monom_basis, monom_function) 
end

export ConstantContinuous

"""
    ConstantContinuous <: AbstractGMPContinuous
   
Type representing constant functions  
"""
struct ConstantContinuous <: AbstractGMPContinuous
    variables::Nothing
    bsa_set::Nothing
    approx_type::EXACT_RELAXATION
    approx_basis::Nothing
    approx_function::Function
end

ConstantContinuous(a::Number) = ConstantContinuous(nothing, nothing, EXACT_RELAXATION(), nothing, x -> x == 1 ? a : 0)
ZeroContinuous() = ConstantContinuous(0)
OneContinuous() = ConstantContinuous(1)

function approximate(f::ConstantContinuous, ::Int)
    return get(f, ApproxFunction())(1)*_mono_one(get(f, Variables()))
end

export VariableContinuous

"""
    VariableContinuous{V, S, R, T}

Type representing a continuous funciton that can be stregthend via a conic approximation. 
"""
struct VariableContinuous{V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet, R <: AbstractApproximation, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::S
    relax_type::R
    relax_basis::T
end
