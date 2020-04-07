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

