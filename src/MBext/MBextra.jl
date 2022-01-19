function MB.change_basis(p::MP.AbstractPolynomialLike, Basis::Type{<:MB.AbstractPolynomialBasis})
    basis = MB.maxdegree_basis(Basis, variables(p), maxdegree(p))
    return MB.change_basis(p, basis)
end

function MB.change_basis(p::MP.AbstractPolynomialLike, basis::MB.AbstractPolynomialBasis)
    coeffs = promote_type(Float64, coefficienttype(p))[]
    rem = p
    mons = monomials(basis)
    for i = 1:length(basis)
        c, rem = divrem(rem, mons[i])
        if iszero(c)
            push!(coeffs, zero(eltype(coeffs)))
        else
            @assert length(monomials(c)) == 1 && MP.isconstant(first(monomials(c)))
            push!(coeffs, first(coefficients(c)))
        end
    end
    idx = findall(!iszero, coeffs)
    return coeffs[idx], mons[idx]
end

function MB.monomials(basis::MB.AbstractPolynomialBasis)
    if basis isa MB.AbstractMonomialBasis
        return basis.monomials
    else
        return basis.polynomials
    end
end
