export MeasureRef
export poly_variables, set_poly_variables
export support, set_support
export relax_type, set_relax_type
export moment_basis, set_moment_basis


struct MeasureRef <: AbstractMeasureRef
    model::GMPModel
    index::Int
end

Base.iszero(::MeasureRef) = false
Base.copy(v::MeasureRef) = MeasureRef(v.model, v.index)
Base.broadcastable(v::MeasureRef) = Ref(v)
JuMP.isequal_canonical(v::MeasureRef, other::MeasureRef) = isequal(v, other)
JuMP.index(v::MeasureRef) = v.index

function Base.hash(v::MeasureRef, h::UInt)
    return hash(objectid(owner_model(v)), hash(index(v), h))
end

function Base.isequal(v1::MeasureRef, v2::MeasureRef)
    return owner_model(v1) === owner_model(v2) && index(v1) == index(v2)
end

function MeasureRef(m::GMPModel)
    index = MOI.add_variable(model_info(m))
    return MeasureRef(m, index)
end

JuMP.name(v::MeasureRef) = model_info(owner_model(v)).variable_names[index(v)]

function JuMP.set_name(v::MeasureRef, s::String)
    model_info(owner_model(v)).variable_names[index(v)] = s
    return v
end

function JuMP.add_variable(model::GMPModel, v::MeasureVariable, name::String="")
    var_ref = MeasureRef(model)
    model_data(model).variables[var_ref] = v
    if !isempty(name)
        set_name(var_ref, name)
    end
    return var_ref   
end

function JuMP.variable_by_name(model::GMPModel, name::String)
    index = findfirst(x -> x==name, model_info(model).variable_names)
    if index isa Nothing
        return nothing
    else
        return MeasureRef(model, index)
    end
end

"""
    _info(v::MeasureRef) 

Get a measure's info attribute.
"""
_info(v::MeasureRef) = model_data(owner_model(v)).variables[v].info

"""
    poly_variables(v::MeasureRef)

Get variables measure is acting on.
"""
poly_variables(v::MeasureRef) = _info(v).poly_variables

"""
    set_poly_variables(v::MeasureRef, pv::Vector{MP.AbstractVariable})

Set variables measure is acting on.
"""
function set_poly_variables(v::MeasureRef, pv::Vector{<:MP.AbstractVariable}) 
    _info(v).poly_variables = pv
    return v
end

"""
    support(v::MeasureRef)

Get support of measure.
"""
support(v::MeasureRef) = _info(v).support

"""
    set_support(v::MeasureRef, set)

Set support of measure.
"""
function set_support(v::MeasureRef, set)
    _info(v).support = set
    return v
end

"""
    moment_basis(v::MeasureRef)

Get basis for moments of measure.
"""
moment_basis(v::MeasureRef) = _info(v).moment_basis

"""
    set_moment_basis(v::MeasureRef, basis)

Set basis for moment of measure.
"""
function set_moment_basis(v::MeasureRef, basis)
    _info(v).moment_basis = basis
    return v
end

"""
    relax_type(v::MeasureRef)

Get relaxation type for measure.
"""
relax_type(v::MeasureRef) = _info(v).relax_type

"""
    set_relax_type(v::MeasureRef, relax_type)

Set relaxation type for measure.
"""
function set_relax_type(v::MeasureRef, relax_type)
    _info(v).relax_type = relax_type
    return v
end
