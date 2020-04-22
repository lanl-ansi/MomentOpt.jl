function SemialgebraicSets.inequalities(::AbstractAlgebraicSet)
    return []
end

function polynomials(s::AbstractSemialgebraicSet)
    return [equalities(s)..., inequalities(s)...]
end

function MultivariatePolynomials.maxdegree(s::AbstractSemialgebraicSet)
    return maximum(maxdegree.(polynomials(s)))
end
