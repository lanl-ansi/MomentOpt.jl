export variable_box, LebesgueMeasure, integrate, measure

# Extend JuMP.VariableRef. Undocumented functions are direct copys of the functions defined for JuMP.VariableRef
abstract type AbstractMeasureRef <: JuMP.AbstractVariableRef end
Base.iszero(::AbstractMeasureRef) = false
Base.copy(v::AbstractMeasureRef) = MeasureRef(v.model, v.index)
Base.broadcastable(v::AbstractMeasureRef) = Ref(v)
JuMP.isequal_canonical(v::AbstractMeasureRef, other::AbstractMeasureRef) = isequal(v, other)
JuMP.owner_model(v::AbstractMeasureRef) = v.model
JuMP.index(v::AbstractMeasureRef) = v.index

function Base.hash(v::AbstractMeasureRef, h::UInt)
    return hash(objectid(owner_model(v)), hash(index(v), h))
end

function Base.isequal(v1::AbstractMeasureRef, v2::AbstractMeasureRef)
    owner_model(v1) === owner_model(v2) && index(v1) == index(v2)
end

JuMP.name(v::AbstractMeasureRef) = model_info(owner_model(v)).variable_names[index(v)]

function JuMP.set_name(v::AbstractMeasureRef, s::String)
    model_info(owner_model(v)).variable_names[index(v)] = s
    return v
end

"""
    measure(v::AbstractMeasureRef)

Get the measure an AbstractMeasureRef is pointing to. 
"""
measure(v::AbstractMeasureRef) = model_data(owner_model(v)).variables[v]

"""
    MP.variables(v::AbstractMeasureRef)

Get variables measure is acting on.
"""
MP.variables(v::AbstractMeasureRef) = MP.variables(measure(v))


"""
    AbstractKnownMeasure


Supertype for all measures whose moments are known explicitly via a function.
For every implementation of a subtype of AbstractKnownMeasures the functions
MP.variables and moment_function have to be defined. 
By default a subtype of AbstractKnownMeasure is assumed to look like
```julia
struct CustomKnownMeasure <: AbstractKnownMeasure
    variables::Vector{MP.AbstractVariable}
    func::Function
end
```
where the variables field is a vector x of variables on which the measure is acting 
and func is a function that takes as input the exponent α of a monomial x^α and retruns
the value of the corresponding moment.
"""
abstract type AbstractKnownMeasure end
MP.variables(m::AbstractKnownMeasure) = m.variables
moment_function(m::AbstractKnownMeasure) = m.func
Base.broadcastable(m::AbstractKnownMeasure) = Ref(m)

"""
    integrate(mon::PT <: MT, meas::AbstractKnownMeasure)

Returns the value of the integral of mon with respect to meas.
"""
function integrate(poly::PT, meas::AbstractKnownMeasure) where PT <: MT 
    vars = MP.variables(meas)
    fun = moment_function(meas)
    if poly isa Number
        c = [poly]
        m = [sum(vars.^0)/length(vars)]
    else
        c = coefficients(poly)
        m = monomials(poly)
    end
    return sum(c[i]*fun(exponents(m[i])) for i in 1:length(c))
end

# In the sequal we show implement the Lebesgue measure on a box. 
# Other measures can be implemented in a similar fasion.

struct VariableBox{T <: Number, VT <: MP.AbstractVariable}
    dict::Dict{VT, Vector{T}}
end
Base.broadcastable(vb::VariableBox) = Ref(vb)
Base.keys(vb::VariableBox) = keys(vb.dict)
"""
    variable_box()

Generates VariableBox that contains information for lower and upper bounds of variables.

## Example

@polyvar x y
vb = variable_box( x=> [0,1], y=>[-1,1])

"""
variable_box(args...) = VariableBox(Dict(args))
 
struct LebesgueMeasure <: AbstractKnownMeasure
    variables::Vector{MP.AbstractVariable}
    func::Function
end

function Base.show(io::IO, m::LebesgueMeasure)
    print(io, "Lebesgue measure acting on: $(MP.variables(m))")
end

lebesgue_line(lb, ub, deg) = (ub^(deg+1) - lb^(deg+1))/(deg+1)

"""
    LebesgueMeasure(variables_with_domain::VariableBox)

Returns an <: AbstractKnonwMeasure representing the Lebesgue measure on the box, defined by the input. 
"""
function LebesgueMeasure(variables_with_domain::VariableBox)
    vars = sort!(collect(keys(variables_with_domain)), rev = true)
    func = e -> prod([lebesgue_line(variables_with_domain.dict[v]..., e[i]) for (i,v) in enumerate(vars)])
    return LebesgueMeasure(vars, func)
end


