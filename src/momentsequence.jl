export moment_sequence, lebesgue_box_sequence, lebesgue_box, variable_box

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
    return MomentSequence(basis, Dict( b => v for (b,v) in zip(collect(basis), values)))
end
MM.moments(ms::MomentSequence) = ms.moments
MP.monomials(ms::MomentSequence) = collect(ms.basis)
MP.maxdegree(ms::MomentSequence) = maximum(maxdegree.(monomials(ms)))
basistype(ms::MomentSequence) = typeof(ms.basis)

function Base.show(io::IO, ms::MomentSequence)
    println(io, "MomentSequence:")
    for (k,v) in  moments(ms)
        print(io, "$k => $v, ")
    end
    println()
end

struct VariableBox{T<:Number, VT<:MP.AbstractVariable}
    dict::Dict{VT, Vector{T}}
end
Base.broadcastable(vb::VariableBox) = Ref(vb)

"""
    variable_box()

Generates VariableBox that contains information for lower and upper bounds of variables.

## Example

@polyvar x y
vb = variable_box( x=> [0,1], y=>[-1,1])

"""
variable_box(args...) = VariableBox(Dict(args))

lebesgue_line(lb, ub, deg) = (ub^(deg+1) - lb^(deg+1))/(deg+1)
function lebesgue_box(
                      variables_with_domain::VariableBox, 
                      f::AbstractMonomial)
    m = 1
    for (v,e) in zip(variables(f), exponents(f))
        m *= lebesgue_line(variables_with_domain.dict[v]..., e)
    end
    return m
end
       
"""
    lebesgue_box(variables_with_domain::VariableBox}, 
                 f::MP.AbstractPolynomialLikei; normalize = false)

Computes the integral of f over the box defined in variables_with_domain. If normalize = true the
lebesgue measure on the box is scaled to be a probability measure.
"""
function lebesgue_box(
                      variables_with_domain::VariableBox, 
                      f::MP.AbstractPolynomialLike; normalize = false)
    mom = sum( c*lebesgue_box(variables_with_domain, m) for (c,m) in zip(coefficients(f), monomials(f)))
    if normalize
        factor = prod( last(r)-first(r) for (k,r) in variables_with_domain.dict)
        return mom/factor
    else
        return mom
    end
end

"""
    lebesgue_box_sequence(
        variables_with_domain::VariableBox,
        max_degree::Int; 
        normalize = false, 
        basis_type = MonomialBasis)

Returns the moment sequence of the lebesgue measure on the box defined by variables_with_domain. 
"""
function lebesgue_box_sequence(variables_with_domain::VariableBox, max_degree::Int; normalize = false, basis_type = MonomialBasis)    
    basis = maxdegree_basis(basis_type, sort!(collect(keys(variables_with_domain.dict)), rev = true), max_degree)
    dict = Dict( mon => lebesgue_box(variables_with_domain, mon; normalize = normalize) for mon in collect(basis))
    return MomentSequence(basis, dict)
end

