function MB.change_basis(p::MP.AbstractPolynomialLike, Basis::Type{<:MB.AbstractPolynomialBasis})
    basis = MB.maxdegree_basis(Basis, variables(p), maxdegree(p))
    return MB.change_basis(p, basis)
end

function MB.change_basis(p::MP.AbstractPolynomialLike, basis::MB.AbstractPolynomialBasis)
    coeffs = Vector{promote_type(Float64, MP.coefficient_type(p))}(undef, length(basis))
    rem = p
    mons = monomials(basis)
    for i in reverse(eachindex(basis))
        c, rem = divrem(rem, mons[i])
        if iszero(c)
            coeffs[i] = zero(eltype(coeffs))
        else
            @assert length(monomials(c)) == 1 && MP.isconstant(first(monomials(c)))
            coeffs[i] = first(coefficients(c))
        end
    end
    idx = findall(!iszero, coeffs)
    return coeffs[idx], mons[idx]
end

function MP.monomials(basis::MB.AbstractPolynomialBasis)
    if basis isa MB.AbstractMonomialBasis
        return basis.monomials
    else
        return basis.polynomials
    end
end
