export approximate

"""
    approximate(f::AbstractGMPContinuous, max_degree::Int)

Returns a polynomial approximation of f of degree max_degree.
"""
function approximate(f::AbstractGMPContinuous, max_degree::Int)
    basis = MB.maxdegree_basis(f, max_degree)
    return dot(eval_vector(basis), basis)
end

export SymbolicContinuous

"""
    SymbolicContinuous{V, S, T} <: AbstractGMPContinuous

Type representing a symbolic continuous function. This type does not allow to compute integrals or relaxations.
"""
struct SymbolicContinuous{S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::S
    approx_type::NO_APPROXIMATION
    approx_basis
end

function SymbolicContinuous(variables::Vector{V}, domain::S, monom_basis) where {S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable} 
    return SymbolicContinuous(support, variables, NO_APPROXIMATION(), monom_basis) 
end

export AnalyticContinuous

"""
    AnalyticContinuous{V, T} <: AbstractGMPContinuous

Type representing an analytic continuous function. Its field coef_function allows to compute the coefficients for arbitrary close polynomial approximations.
"""
struct AnalyticContinuous{V <: MP.AbstractVariable} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::Nothing
    approx_type::EXACT_APPROXIMATION
    approx_basis
    approx_function::Function
end

function AnalyticContinuous(variables::Vector{V}, monom_basis, monom_function::Function) where {V <: MP.AbstractVariable} 
    return AnalyticContinuous(variables, nothing,  EXACT_APPROXIMATION(), monom_basis, monom_function) 
end

export ConstantContinuous, ZeroContinuous, OneContinuous

"""
    ConstantContinuous <: AbstractGMPContinuous
   
Type representing constant functions  
"""
struct ConstantContinuous <: AbstractGMPContinuous
    variables::Nothing
    bsa_set::Nothing
    approx_type::EXACT_APPROXIMATION
    approx_basis::Nothing
    approx_function::Function
end

ConstantContinuous(a::Number) = ConstantContinuous(nothing, nothing, EXACT_APPROXIMATION(), nothing, x -> x == 1 ? a : 0)
ZeroContinuous() = ConstantContinuous(0)
OneContinuous() = ConstantContinuous(1)

function approximate(f::ConstantContinuous, ::Int)
    return MOI.get(f, ApproxFunction())(1)*_mono_one(get(f, Variables()))
end

export VariableContinuous

"""
    VariableContinuous{V, S, R}

Type representing a continuous funciton that can be stregthend via a conic approximation. 
"""
struct VariableContinuous{V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet, R <: AbstractApproximation} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::S
    approx_type::R
    approx_basis
end


export Cont
Cont(vars; domain = FullSpace(), basis = MonomialBasis, approx = DefaultApproximation()) = VariableContinuous(vars, domain, approx, basis)

