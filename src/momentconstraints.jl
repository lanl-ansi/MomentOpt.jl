export MomentConstraint

struct MomentConstraintShape <: JuMP.AbstractShape end
Base.broadcastable(shape::MomentConstraintShape) = Ref(shape)

"""
    MomentConstraint

Type representing a moment constraint, i.e., a constraint setting a moment expression in relation to a number.     
"""
mutable struct MomentConstraint{T <: Number, PT<:MT} <: AbstractGMPConstraint
	func::MomentExpr{T, PT}
    set::MOI.AbstractScalarSet
end

function MomentConstraint(func::AffineMomentExpr, set::Union{MOI.EqualTo{T}, MOI.LessThan{T}, MOI.GreaterThan{T}}) where T
    return MomentConstraint(momexpr(func), MOI.Utilities.shift_constant(set, convert(T,-constant(func))))
end

function Base.promote_rule(::Type{MomentConstraint{T1, PT1}},::Type{MomentConstraint{T2, PT2}})  where {T1 <: Number, T2 <: Number, PT1<:MT, PT2<:MT}
    return MomentConstraint{promote_type(T1, T2), promote_type(PT1, PT2)}
end

function Base.convert(::Type{MomentConstraint{T, PT1}}, mc::MomentConstraint) where {T <: Number, PT1 <: MT}
    return MomentConstraint{T, PT1}(convert(MomentExpr{T, PT1}, mc.func), mc.set)
end

JuMP.jump_function(con::MomentConstraint) = con.func
JuMP.moi_set(con::MomentConstraint) = con.set
JuMP.shape(con::MomentConstraint) = MomentConstraintShape()
JuMP.reshape_vector(expr::MomentExpr, ::MomentConstraintShape) = expr
JuMP.reshape_set(set::MOI.AbstractScalarSet, ::MomentConstraintShape) = set

function Base.show(io::IO, mc::MomentConstraint)
    print(io, JuMP.constraint_string(JuMP.REPLMode, mc))
end

function JuMP.function_string(mode, mc::MomentConstraint)
    return string(mc.func)
end

function measures(mc::MomentConstraint)
    return measures(mc.func)
end

function measures(mcv::Array{MomentConstraint{T}}) where T<: MT
    measv = Set{Measure}()
	for mc in mcv
        union!(measv, measures(mc))
	end
    return collect(measv)
end

function constant(mc::MomentConstraint)
    return MOI.constant(mc.set)
end
