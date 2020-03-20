export MeasureRef
export poly_variables, set_poly_variables
export support, set_support
export relax_type, set_relax_type
export moment_basis, set_moment_basis


struct MeasureRef <: AbstractMeasureRef
    model::GMPModel
    index::Int
end

function MeasureRef(m::GMPModel)
    index = MOI.add_variable(model_info(m))
    return MeasureRef(m, index)
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
    index = findfirst(x -> x == name, model_info(model).variable_names)
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
_info(v::MeasureRef) = info(measure(v))

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

struct KnownMeasureRef <: AbstractMeasureRef
    model::GMPModel
    index::Int
end

function KnownMeasureRef(m::GMPModel)
    index = MOI.add_variable(model_info(m))
    return KnownMeasureRef(m, index)
end

function JuMP.add_variable(model::GMPModel, v::KnowMeasureVariable, name::String="")
    var_ref = MeasureRef(model)
    model_data(model).known_measures[var_ref] = v
    if !isempty(name)
        set_name(var_ref, name)
    end
    return var_ref   
end

function known_measure_by_name(model::GMPModel, name::String)
    index = findfirst(x -> x == name, model_info(model).variable_names)
    if index isa Nothing
        return nothing
    else
        return KnownMeasureRef(model, index)
    end
end

"""
    integrate(mon::PT <: MT, meas::KnownMeasureRef)

Returns the value of the integral of mon with respect to meas.
"""
function integrate(poly::PT, meas::KnownMeasureRef) where PT <: MT 
    return integrate(poly, measure_variable(meas))
end
