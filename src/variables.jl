abstract type AbstractGMPVariable  <: JuMP.AbstractVariable end
abstract type AbstractGMPVariableRef <: JuMP.AbstractVariableRef end

export GMPVariable

"""
    GMPVariable <: AbstractGMPVariable

Type representing a variable for the GMPModel.
"""
struct GMPVariable <: AbstractGMPVariable
    v::VariableMeasure
end

object(v::GMPVariable) = v.v
object_type(v) = typeof(object(v))

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, m::VariableGMPObject; extra_kwargs...)
    return GMPVariable(m)
end

# define GMPVariableRef

export GMPVariableRef
"""
    GMPVariableRef

Links a variable to a model.
"""
struct GMPVariableRef <: AbstractGMPVariableRef
    model::JuMP.AbstractModel
    index::Int
end

# define internal functions for GMPVariableRef
Base.broadcastable(v::GMPVariableRef) = Ref(v)
Base.iszero(::GMPVariableRef) = false
JuMP.isequal_canonical(v::GMPVariableRef, w::GMPVariableRef) = v == w
JuMP.index(vref::GMPVariableRef) = vref.index
Base.:(==)(v::GMPVariableRef, w::GMPVariableRef) = v.model === w.model && v.index == w.index
