function SemialgebraicSets.inequalities(::AbstractAlgebraicSet)
    return []
end

function polynomials(s::AbstractSemialgebraicSet)
    return [equalities(s)..., inequalities(s)...]
end

function MultivariatePolynomials.maxdegree(s::AbstractSemialgebraicSet)
    pols = polynomials(s)
    return isempty(pols) ? 0 : maximum(maxdegree.(pols))
end
