export GMPVariable

"""
    GMPVariable{T <: AbstractGMPObject} <: AbstractGMPVariable

Type representing a variable for the GMPModel.
"""
struct GMPVariable{S, T} <: AbstractGMPVariable
    v::S
    function GMPVariable(o::AbstractGMPObject)
        S = typeof(o)
        T = supertype(S)
        new{S,T}(o)
    end
end

object_type(::GMPVariable{S, T}) where {S, T} = S
variable_type(::GMPVariable{S, T}) where {S, T} = T

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, m::AbstractGMPObject; extra_kwargs...)
    return GMPVariable(m)
end

# define GMPVariableRef

"""
    GMPVariableRef

Links a variable to a model.
"""
struct GMPVariableRef <: AbstractGMPVariableRef
    model::JuMP.AbstractModel
    index::Int
    var_type
end

# define internal functions for GMPVariableRef
Base.broadcastable(v::GMPVariableRef) = Ref(v)
Base.iszero(::GMPVariableRef) = false
JuMP.isequal_canonical(v::GMPVariableRef, w::GMPVariableRef) = v == w
Base.:(==)(v::GMPVariableRef, w::GMPVariableRef) = v.model === w.model && v.index == w.index
Base.copy(v::GMPVariableRef) = v
Base.copy(v::GMPVariableRef, m::JuMP.AbstractModel) = GMPVariableRef(m, v.index)
variable_ref_type(vref::GMPVariableRef) = vref.var_type

"""
    compatible(vref1::GMPVariableRef, vref2::GMPVariableRef) 

Checks whether two vrefs can be combined in the same GMPAffexpr.
"""
compatible(vref1::GMPVariableRef, vref2::GMPVariableRef) = variable_ref_type(vref1) == variable_ref_type(vref2)
