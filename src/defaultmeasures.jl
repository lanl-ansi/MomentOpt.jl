export ZeroMeasure

"""
    ZeroMeasure(vars::Vector{MP.AbstractVariable})

Measure acting on vars and giving zero mass to the whole space. 
"""
ZeroMeasure(vars::Vector{<:MP.AbstractVariable}; basis = MonomialBasis) = AnalyticMeasure(vars, basis, x -> 0)

export DiracMeasure

"""
    DiracMeasure(vars::Vector{MP.AbstractVariable}, point::AbstractVector)

Dirac measure in the point vars == point.
"""
function DiracMeasure(vars::Vector{<:MP.AbstractVariable}, point::AbstractVector)
    @assert length(vars) == length(point) "Variables and point need to have same length."
    return AnalyticMeasure(vars, MonomialBasis, x -> prod(point[i]^exponents(x)[i] for i in 1:length(point)))
end

struct VariableBox{T<:Number, VT<:MP.AbstractVariable}
    dict::Dict{VT, Vector{T}}
end
Base.broadcastable(vb::VariableBox) = Ref(vb)

export variable_box

"""
    variable_box()

Generates VariableBox that contains information for lower and upper bounds of variables.

## Example

@polyvar x y
vb = variable_box( x=> [0,1], y=>[-1,1])

"""
function variable_box(args...)
    T = Int
    for i in 1:length(args)
        T = promote_type(T, eltype(last(args[i])))
    end
    return VariableBox(Dict(first(a) => convert(Vector{T}, last(a)) for a in args))
end

lebesgue_line(lb, ub, deg) = (ub^(deg+1) - lb^(deg+1))/(deg+1)
function lebesgue_box(
                      variables_with_domain::VariableBox, 
                      f::AbstractMonomial)
    m = 1
    for (v, e) in zip(variables(f), exponents(f))
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
    mom = sum( c*lebesgue_box(variables_with_domain, m) for (c, m) in zip(coefficients(f), monomials(f)))
    if normalize
        factor = prod( last(r)-first(r) for (k, r) in variables_with_domain.dict)
        return mom/factor
    else
        return mom
    end
end

export lebesgue_measure_box
"""
    lebesgue_measure_box(
        variables_with_domain::VariableBox,
        normalize = false, 
        basis = MonomialBasis)

Returns the Lebesgue measure on the box defined by variables_with_domain. 
"""
function lebesgue_measure_box(variables_with_domain::VariableBox; normalize = false, basis = MonomialBasis)    
    return AnalyticMeasure(
                           sort!(collect(keys(variables_with_domain.dict)), rev = true),
                           basis, 
                           x -> lebesgue_box(variables_with_domain, x; normalize = normalize)
                          )
end


# TODO add gaussian
