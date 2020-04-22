function SemialgebraicSets.inequalities(::AbstractAlgebraicSet)
    return []
end

function polynomials(s::AbstractSemialgebraicSet)
    return [equalities(s)..., inequalities(s)...]
end

function MultivariatePolynomials.maxdegree(s::AbstractSemialgebraicSet)
    pol = polynomials(s)
    return isempty(pol) ? 0 : maximum(maxdegree.(pols))
end
