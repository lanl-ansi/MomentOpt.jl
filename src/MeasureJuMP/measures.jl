export integrate

"""
     integrate(p::Number, m::AbstractGMPMeasure)
     integrate(p::MP.AbstractPolynomialLike, m::AbstractGMPMeasure)

Returns the integral of p with respect to m. 
"""
function integrate(p::Number, m::AbstractGMPMeasure)
    return p*eval_vector(MB.maxdegree_basis(m, 0))
end

function integrate(p::MP.AbstractPolynomialLike, m::AbstractGMPMeasure)
    basis = covering_basis(m, p)
    coefs, _ = MB.change_basis(p, basis)
    return dot(coefs, eval_vector(basis, m))
end

export SymbolicMeasure

"""
    SymbolicMeasure{V,S, T} <: AbstractGMPMeasure

Type representing a symbolic measure. This type does not allow to compute integrals or relaxations.
"""
struct SymbolicMeasure{ V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet} <: AbstractGMPMeasure
    variables::Vector{V}
    bsa_set::S
    approx_type::NO_APPROXIMATION
    approx_basis
end

function SymbolicMeasure(variables::Vector{V}, support::S, moment_basis::T) where { V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet, T <: Type{<: MB.AbstractPolynomialBasis}}
    return SymbolicMeasure(variables, support, NO_APPROXIMATION(), moment_basis) 
end

export AnalyticMeasure

"""
    AnalyticMeasure{V, T} <: AbstractGMPMeasure

Type representing an analytic measure. Its field approx_function allows to integrate polynomials."""
struct AnalyticMeasure{V <: MP.AbstractVariable} <: AbstractGMPMeasure
    variables::Vector{V}
    bsa_set::Nothing
    approx_type::EXACT_APPROXIMATION
    approx_basis
    approx_function::Function
end

function AnalyticMeasure(variables::Vector{V}, moment_basis::T, moment_function::Function) where  {V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} 
    return AnalyticMeasure(variables, nothing, EXACT_APPROXIMATION(), moment_basis, moment_function)
end

export ZeroMeasure

"""
    ZeroMeasure <: AbstractGMPMeasure
   
Type representing the measure giving zero weight to anything. 
"""
struct ZeroMeasure <: AbstractGMPMeasure
    variables::Nothing
    bsa_set::Nothing
    approx_type::EXACT_APPROXIMATION
    approx_basis::Nothing
    approx_function::Function
end

ZeroMeasure() = ZeroMeasure(nothing, nothing, EXACT_APPROXIMATION(), nothing, x -> 0)
Base.show(io::IO, ::ZeroMeasure) = show(io, 0) 

function covering_basis(t::ZeroMeasure, p::MP.AbstractPolynomialLike)
    return MB.basis_covering_monomials(MonomialBasis, monomials(p)) 
end

function integrate(p::Number, ::ZeroMeasure)
    return zero(typeof(p))
end

function integrate(p::MP.AbstractPolynomialLike, ::ZeroMeasure)
    return zero(coefficienttype(p))
end

export VariableMeasure

"""
    VariableMeasure{V, S, R, T}

Type representing a measure that can be relaxed via a conic relaxation. 
"""
struct VariableMeasure{V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet, R <: AbstractApproximation} <: AbstractGMPMeasure
    variables::Vector{V}
    bsa_set::S
    approx_type::R
    approx_basis
end

export Meas
Meas(vars; support = FullSpace(), basis = MonomialBasis, approx = DefaultApproximation()) = VariableMeasure(vars, support, approx, basis)

