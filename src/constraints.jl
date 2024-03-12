abstract type AbstractGMPSet <: MOI.AbstractScalarSet end
abstract type AbstractGMPShape <: JuMP.AbstractShape end
abstract type AbstractGMPConstraint <: JuMP.AbstractConstraint end

function Base.:(==)(c1::AbstractGMPConstraint, c2::AbstractGMPConstraint)
    return JuMP.jump_function(c1) == JuMP.jump_function(c2) && JuMP.moi_set(c1) == JuMP.moi_set(c2)
end

# MomentConstraint
"""
    MomentConstraintShape 

"""
struct MomentConstraintShape <: AbstractGMPShape end

JuMP.reshape_vector(expr::MomentExpr, ::MomentConstraintShape) = expr
JuMP.reshape_set(set::MOI.AbstractScalarSet, ::MomentConstraintShape) = set
#TODO dual shape for MomentConstraintShape

"""
    MomentConstraint

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

JuMP.in_set_string(mode, m::AnalyticMeasure) = "= "*sprint(show, m)

function Base.show(io::IO, con::MomentConstraint)
    print(io, JuMP.constraint_string(MIME"text/plain"(), con))
end

"""
    MeasureConstraintShape 

"""
struct MeasureConstraintShape <: AbstractGMPShape end

JuMP.reshape_vector(expr::MeasExpr, ::MeasureConstraintShape) = expr
JuMP.reshape_set(set::AnalyticMeasure, ::MeasureConstraintShape) = set

#TODO dual shape for MeasureConstraintShape


"""
    MeasureConstraint

"""
struct MeasureConstraint <: AbstractGMPConstraint
    func::MeasExpr
    set::AnalyticMeasure
end

MeasureConstraint(o::MeasExpr) = MeasureConstraint(o, ZeroMeasure(variables(o)))
MeasureConstraint(ae::AffMeasExpr) = MeasureConstraint(expr(ae), -constant(ae))

JuMP.jump_function(con::MeasureConstraint) = con.func
JuMP.moi_set(con::MeasureConstraint) = con.set
JuMP.shape(con::MeasureConstraint) = MeasureConstraintShape()

function JuMP.function_string(mode, mc::MeasureConstraint)
    return string(mc.func)
end

function Base.show(io::IO, con::MeasureConstraint)
    print(io, JuMP.constraint_string(MIME"text/plain"(), con))
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
#Base.copy(v::GMPConstraintRef) = v
#Base.copy(v::GMPConstraintRef, m::JuMP.AbstractModel) = GMPConstraintRef(m, v.index, v.shape)
JuMP.owner_model(cref::GMPConstraintRef) = cref.model
JuMP.shape(cref::GMPConstraintRef) = cref.shape
JuMP.index(cref::GMPConstraintRef) = cref.index
cref_object(cref::GMPConstraintRef) = gmp_constraints(owner_model(cref))[index(cref)]
JuMP.dual(vref::GMPConstraintRef) = JuMP.dual(vref, shape(vref))
