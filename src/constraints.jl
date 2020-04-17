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

function Base.show(io::IO, con::MomentConstraint)
    print(io, JuMP.constraint_string(JuMP.REPLMode, con))
end

function JuMP.build_constraint(_error::Function, ae::AffMomentExpr, set::MOI.AbstractScalarSet)
    return MomentConstraint(ae, set)
end

# MeasureConstraints
"""
    EqualToMeasure <: AbstractGMPSet

The set consiting of only a single measure.
"""
struct EqualToMeasure <: AbstractGMPSet
    measure::AnalyticMeasure
end

MOI.constant(set::EqualToMeasure) = set.measure
MOI.Utilities.shift_constant(set::EqualToMeasure, m::AnalyticMeasure) = EqualToMeasure(constant(set) + m)
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


"""
    MeasureConstraint

"""
struct MeasureConstraint <: AbstractGMPConstraint
    func::ObjectExpr
    set::EqualToMeasure
        function MeasureConstraint(func::ObjectExpr, set::EqualToMeasure)
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
    return MeasureConstraint(expr(ae), EqualToMeasure(constant(ae))) 
end


function JuMP.build_constraint(_error::Function, e::ObjectExpr, s::EqualToMeasure)
    return MeasureConstraint(e, s) 
end

function JuMP.build_constraint(_error::Function, ae::AffObjectExpr, s::EqualToMeasure)
    return MeasureConstraint(expr(ae), MOIU.shift_constant(s, -constant(ae))) 
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
