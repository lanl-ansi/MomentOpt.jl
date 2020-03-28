export GMPVariable

"""
    GMPVariable{T <: AbstractGMPObject} <: AbstractGMPVariable

Type representing a variable for the GMPModel.
"""
struct GMPVariable{T <: AbstractGMPObject} <: AbstractGMPVariable
    v::T
end

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, m::AbstractGMPObject; extra_kwargs...)
    return GMPVariable(m)
end

# define GMPVariableRef

"""
    GMPVariableRef

Links a variable to a model.
"""
struct GMPVariableRef{T <: JuMP.AbstractModel, V } <: AbstractGMPVariableRef
    model::T
    index::MOI.VariableIndex
    var_type::V
end

# define internal functions for GMPVariableRef
Base.broadcastable(v::GMPVariableRef) = Ref(v)
Base.iszero(::GMPVariableRef) = false
JuMP.isequal_canonical(v::GMPVariableRef, w::GMPVariableRef) = v == w
Base.:(==)(v::GMPVariableRef, w::GMPVariableRef) = v.model === w.model && v.index == w.index
Base.copy(v::GMPVariableRef) = v
Base.copy(v::GMPVariableRef, m::JuMP.AbstractModel) = GMPVariableRef(m, v.index)

"""
    compatible(vref1::GMPVariableRef, vref2::GMPVariableRef) 

Checks whether two vrefs can be combined in the same GMPAffexpr.
"""
compatible(vref1::GMPVariableRef{S}, vref2::GMPVariableRef{T}) = S == T

