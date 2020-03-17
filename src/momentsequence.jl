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
MP.monomials(ms::MomentSequence) = collect(ms.basis)
MP.maxdegree(ms::MomentSequence) = maximum(maxdegree.(monomials(ms)))
basistype(ms::MomentSequence) = typeof(ms.basis)

function Base.show(io::IO, ms::MomentSequence)
    print(io, moments(ms))
end

lebesgue_line(lb, ub, deg) = (ub^(deg+1) - lb^(deg+1))/(deg+1)
function lebesgue_box(
                      variables_with_domain::Dict{VT, Vector{T}}, 
                      f::Monomial) where {VT<:MP.AbstractVariable, T<:Number}
    m = 1
    for (v,e) in zip(variables(f), exponents(f))
        m *= lebesgue_line(variables_with_domain[v]..., e)
    end
    return m
end
       
"""
    lebesgue_box(variables_with_domain::Dict{MP.AbstractVariable, Vector{<:Number}}, 
                 f::MP.AbstractPolynomialLikei; normalize = false)

Computes the integral of f over the box defined in variables_with_domain. If normalize = true the
lebesgue measure on the box is scaled to be a probability measure.
"""
function lebesgue_box(
                      variables_with_domain::Dict{VT, Vector{T}}, 
                      f::MP.AbstractPolynomialLike; normalize = false) where {VT<:MP.AbstractVariable, T<:Number}
    mom = sum( c*lebesgue_box(variables_with_domain, m) for (c,m) in zip(coefficients(f), monomials(f)))
    if normalize
        factor = prod( last(r)-first(r) for (k,r) in variables_with_domain)
        return mom/factor
    else
        return mom
    end
end

"""
    lebesgue_box_sequence(
        variables_with_domain::Dict{VT, Vector{T}},
        max_degree::Int; 
        normalize = false, 
        basis_type = MonomialBasis) where {VT<:MP.AbstractVariable, T<:Number} 

Returns the moment sequence of the lebesgue measure on the box defined by variables_with_domain. 
"""
function lebesgue_box_sequence(variables_with_domain::Dict{VT, Vector{T}},max_degree::Int; normalize = false, basis_type = MonomialBasis) where {VT<:MP.AbstractVariable, T<:Number}    
    basis = maxdegree_basis(basis_type, sort!(collect(keys(variables_with_domain)), rev = true), max_degree)
    dict = Dict( mon => lebesgue_box(variables_with_domain, mon; normalize = normalize) for mon in collect(basis))
    return MomentSequence(basis, dict)
end

