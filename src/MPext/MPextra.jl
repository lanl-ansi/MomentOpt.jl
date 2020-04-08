function MP.monomials(basis::MB.AbstractPolynomialBasis)
    if basis isa MB.AbstractMonomialBasis
        return basis.monomials
    else
        return basis.polynomials
    end
end

MP.maxdegree(::Number) = 0
