using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases
using MultivariateMoments
const MM = MultivariateMoments

using DynamicPolynomials

Base.collect(basis::MB.AbstractPolynomialVectorBasis) = basis.polynomials
Base.collect(basis::MB.AbstractMonomialBasis) = basis.monomials
MP.monomials(
             basis_type::Type{<:MB.AbstractPolynomialBasis}, 
             vars::AbstractVector{PT}, 
             range::UnitRange{Int}) where PT <: MP.AbstractVariable = 
collect(maxdegree_basis(basis_type, vars, last(range)))

"""
    BasisPolynomial{T<:Number, BT <: MB.AbstractPolynomialBasis}

Represends polynomial in a particular basis.
"""
mutable struct BasisPolynomial{T, BT <: MB.AbstractPolynomialBasis}
    coeffs::Vector{T}
    basis::BT
end

MP.coefficients(f::BasisPolynomial) = f.coeffs
MP.monomials(f::BasisPolynomial) = MP.monomials(f.basis)

function Base.show(io::IO, f::BasisPolynomial)
    l = length(coefficients(f))
    for i = 1:l-1
        print(io, "$(coefficients(f)[i]) * ( $(monomials(f)[i]) ) + ")
    end
    print(io, "$(coefficients(f)[l]) * ( $(monomials(f)[l]) )")
end

"""
    basis_polynomial(basis_type::Type{<:MB.AbstractPolynomialBasis}, f::MP.AbstractPolynomialLike)

Converts a MultivariatePolynomials polynomial f into a BasisPolynomial in basis_type.
"""
function basis_polynomial(basis_type::Type{<:MB.AbstractPolynomialBasis}, f::MP.AbstractPolynomialLike)
    r = f
    coeffs = []
    basis = basis_covering_monomials(basis_type, monomials(f))
    for p in collect(basis)
        c, r = divrem(r, p)
        push!(coeffs, c)
    end
    return BasisPolynomial(coeffs, basis)
end

MP.polynomial(f::BasisPolynomial) = MP.polynomial(i -> coefficients(f)[i], f.basis)

"""
    MomentSequence{T, PT <: MP.AbstractPolynomialLike, BT<:Type{<:MB.AbstractPolynomialBasis}}

Represends a truncated moment sequence in a particular basis.
"""
struct MomentSequence{T, PT <: MP.AbstractPolynomialLike, BT<:MB.AbstractPolynomialBasis}
    basis::BT
    moments::Dict{PT, T}
end

Base.broadcastable(ms::MomentSequence) = Ref(ms)
"""
    moment_sequence(basis::MB.AbstractPolynomialBasis, values::Vector{T})
Constructs a MomentSequence. 
"""
function moment_sequence(basis::MB.AbstractPolynomialBasis, values::Vector{T}) where T
    @assert length(basis) == length(values)
    return MomentSequence(basis, Dict( b => v for (b,v) in zip(monomials(basis), values)))
end
MM.moments(ms::MomentSequence) = ms.moments
MP.monomials(ms::MomentSequence) = monomials(ms.basis)
MP.maxdegree(ms::MomentSequence) = maximum(maxdegree.(monomials(ms)))
basistype(ms::MomentSequence) = typeof(ms.basis)

function Base.show(io::IO, ms::MomentSequence)
    print(io, moments(ms))
end

"""
    riesz(ms::MomentSequence, f::MP.AbstractPolynomialLike)

Applies the Riesz functional associated to ms to f.
"""
function riesz(ms::MomentSequence, f::MP.AbstractPolynomialLike)
    @assert maxdegree(f) <= maxdegree(ms)
    fb = basis_polynomial(basistype(ms), f)
    @assert monomials(fb) ⊆ monomials(ms)
    return sum( coef*moments(ms)[mon] for (coef, mon) in zip(coefficients(fb), monomials(fb)))
end


