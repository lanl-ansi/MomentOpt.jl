function MP.monomials(basis::MB.AbstractPolynomialBasis)
    if basis isa MB.AbstractMonomialBasis
        return basis.monomials
    else
        return basis.polynomials
    end
end

function MB.change_basis(p::MP.AbstractPolynomialLike, Basis::Type{<:MB.AbstractPolynomialBasis})
    basis = MB.maxdegree_basis(Basis, variables(p), maxdegree(p))
    coeffs = Float64[]
    rem = p
    mons = monomials(basis)
    for i = 1:length(basis)
        c, rem = divrem(rem, mons[i])
        push!(coeffs, c)
    end
    idx = findall(!iszero, coeffs)
    return coeffs[idx], mons[idx]
end

function MB.change_basis(p::MP.AbstractPolynomialLike, basis::MB.AbstractPolynomialBasis)
    coeffs = Float64[]
    rem = p
    mons = monomials(basis)
    for i = 1:length(basis)
        c, rem = divrem(rem, mons[i])
        push!(coeffs, c)
    end
    idx = findall(!iszero, coeffs)
    return coeffs[idx], mons[idx]
end

MP.maxdegree(::Number) = 0
