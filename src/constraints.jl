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

function Base.show(io::IO, con::MomentConstraint)
    print(io, JuMP.constraint_string(JuMP.REPLMode, con))
end

function JuMP.build_constraint(_error::Function, ae::AffMomentExpr, set::MOI.AbstractScalarSet)
    return MomentConstraint(ae, set)
end

# MeasureConstraints
"""
    EqualToMeasure{T} <: AbstractScalarSet

The set consiting of only a single measure.
"""
struct EqualToMeasure{T <: AbstractGMPMeasure} <: AbstractGMPSet
    measure::T
end

MOI.constant(set::EqualToMeasure) = set.measure

function JuMP.in_set_string(print_mode, set::EqualToMeasure)
    return "= "*sprint(show, MOI.constant(set))
end

# TODO MOI.dual_set for EqualToMeasure (GreaterThenContinuous)

"""
    MeasureConstraintShape 

"""
struct MeasureConstraintShape <: AbstractGMPShape end

JuMP.reshape_vector(expr::ObjectExpr, ::MeasureConstraintShape) = expr
JuMP.reshape_set(set::EqualToMeasure, ::MeasureConstraintShape) = set
#TODO dual shape for MeasureConstraintShape

function same_basis(mv::Vector{GMPVariableRef})
    if isempty(mv)
        return
    else
        basis = approx_basis(first(mv))
        if all( x -> x == basis, approx_basis.(mv[2:end]))
            return
        else
            @warn("Measures have different ApproximationBasis. $basis will be used for MeasureConstraint.")
            return 
        end
    end
end
"""
    MeasureConstraint

"""
struct MeasureConstraint <: AbstractGMPConstraint
    func::ObjectExpr
    set::EqualToMeasure
        function MeasureConstraint(func::ObjectExpr, set::EqualToMeasure)
            same_basis(gmp_variables(func))
            new(func, set)
        end
end

JuMP.jump_function(con::MeasureConstraint) = con.func
JuMP.moi_set(con::MeasureConstraint) = con.set
JuMP.shape(con::MeasureConstraint) = MeasureConstraintShape()


function JuMP.function_string(mode, mc::MeasureConstraint)
    return string(mc.func)
end

function Base.show(io::IO, con::MeasureConstraint)
    print(io, JuMP.constraint_string(JuMP.REPLMode, con))
end

function JuMP.build_constraint(_error::Function, ae::AffObjectExpr, s::MOI.EqualTo)
    @assert(constant(ae) == 0, "Affine measure constraints are not supported.")
    @assert(MOI.constant(s) == 0, "Affine measure constraints are not supported.")
    return MeasureConstraint(expr(ae), EqualToMeasure(ZeroMeasure())) 
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

