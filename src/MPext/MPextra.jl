MP.maxdegree(::Number) = 0

"""
    polypoly(p::AbstractPolynomialLike, vars::Vector{MP.AbstractVariable})

Retruns a pair of vectors (coefs, mons), where coefs is a vector of AbstractPolynomialLike and mons is a vector of monomials in vars such that dot(coef, mons) == p 
"""
function polypoly(p::AbstractPolynomialLike, vars::Vector{<:MP.AbstractVariable})
    @assert vars ⊆ variables(p) 
    coefs = coefficients(p).*[subs(monomials(p)[i], [v => 1 for v in vars]...) for i=1:length(monomials(p))]
    mons = subs(monomials(p), [v => 1 for v in setdiff(variables(p), vars)]...)
    return coefs, mons
end
