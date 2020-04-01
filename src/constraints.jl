# MomentConstraint
"""
    MomentConstraintShape 

"""
struct MomentConstraintShape <: AbstractGMPShape end

JuMP.reshape_vector(expr::MomentExpr, ::MomentConstraintShape) = expr
JuMP.reshape_set(set::MOI.AbstractScalarSet, ::MomentConstraintShape) = set
#TODO dual shape for MomentConstraintShape

"""
    MomentConstraint{}

"""
struct MomentConstraint <: AbstractGMPConstraint
    func::MomentExpr
    set::MOI.AbstractScalarSet
end

function MomentConstraint(ae::AffMomentExpr, set)
    return MomentConstraint(ae.expr, MOI.Utilities.shift_constant(set, -convert(typeof(MOI.constant(set)), ae.cons)))
end

JuMP.jump_function(con::MomentConstraint) = con.func
JuMP.moi_set(con::MomentConstraint) = con.set
JuMP.shape(con::MomentConstraint) = MomentConstraintShape()

function JuMP.function_string(mode, mc::MomentConstraint)
    return string(mc.func)
end

function Base.show(io::IO, con::MomentConstraint)
    print(io, JuMP.constraint_string(JuMP.REPLMode, con))
end

function JuMP.build_constraint(_error::Function, ae::AffMomentExpr, set::MOI.AbstractScalarSet)
    return MomentConstraint(ae, set)
end

"""
    GMPConstraintRef.

Links a constraint to a model
"""
abstract type AbstractGMPConstraintRef end

struct GMPConstraintRef <: AbstractGMPConstraintRef
    model::JuMP.AbstractModel
    index::Int
    shape::AbstractGMPShape
end

Base.broadcastable(v::GMPConstraintRef) = Ref(v)
Base.iszero(::GMPConstraintRef) = false
JuMP.isequal_canonical(v::GMPConstraintRef, w::GMPConstraintRef) = v == w
Base.:(==)(v::GMPConstraintRef, w::GMPConstraintRef) = v.model === w.model && v.index == w.index
Base.copy(v::GMPConstraintRef) = v
Base.copy(v::GMPConstraintRef, m::JuMP.AbstractModel) = GMPConstraintRef(m, v.index, v.shape)
JuMP.owner_model(con_ref::GMPConstraintRef) = con_ref.model
JuMP.shape(cref::GMPConstraintRef) = cref.shape
JuMP.index(cref::GMPConstraintRef) = cref.index

