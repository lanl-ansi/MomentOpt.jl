abstract type AbstractGMPVariable  <: JuMP.AbstractVariable end
abstract type AbstractGMPVariableRef <: JuMP.AbstractVariableRef end

export GMPVariable

"""
    GMPVariable{S, T <: AbstractGMPObject} <: AbstractGMPVariable

Type representing a variable for the GMPModel.
"""
struct GMPVariable{S} <: AbstractGMPVariable
    v::S
    function GMPVariable(o::VariableGMPObject)
        S = supertype(typeof(o))
        new{S}(o)
    end
end

variable_object(v::GMPVariable) = v.v
object_type(v::GMPVariable) = typeof(variable_object(v))

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, m::VariableGMPObject; extra_kwargs...)
    return GMPVariable(m)
end

# define GMPVariableRef

export GMPVariableRef
"""
    GMPVariableRef

Links a variable to a model.
"""
struct GMPVariableRef{S <: AbstractGMPObject} <: AbstractGMPVariableRef
    model::JuMP.AbstractModel
    index::Int
    var_type::Type{S}
end

# define internal functions for GMPVariableRef
Base.broadcastable(v::GMPVariableRef) = Ref(v)
Base.iszero(::GMPVariableRef) = false
JuMP.isequal_canonical(v::GMPVariableRef, w::GMPVariableRef) = v == w
JuMP.index(vref::GMPVariableRef) = vref.index
vref_type(vref::GMPVariableRef{S}) where S = S
Base.:(==)(v::GMPVariableRef, w::GMPVariableRef) = v.model === w.model && v.index == w.index && vref_type(v) == vref_type(w)
