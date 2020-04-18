function MB.change_basis(p::MP.AbstractPolynomialLike, Basis::Type{<:MB.AbstractPolynomialBasis})
    basis = MB.maxdegree_basis(Basis, variables(p), maxdegree(p))
    return change_basis(p, basis)
end

function MB.change_basis(p::MP.AbstractPolynomialLike, basis::MB.AbstractPolynomialBasis)
    coeffs = promote_type(Float64, coefficienttype(p))[]
    rem = p
    mons = monomials(basis)
    for i = 1:length(basis)
        c, rem = divrem(rem, mons[i])
        if c == 0
            push!(coeffs, 0)
        else
            @assert length(monomials(c)) == 1 && MP.isconstant(first(monomials(c)))
            push!(coeffs, first(coefficients(c)))
        end
    end
    idx = findall(!iszero, coeffs)
    return coeffs[idx], mons[idx]
end
