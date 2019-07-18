export MomObj, MomCon, MomCons
# export NLMomObj

abstract type AbstractMomentObjective end
struct EmptyObjective <: AbstractMomentObjective end
struct MomConShape <: JuMP.AbstractShape end

"""
MomObj
"""
mutable struct MomObj{PT<:MT} <:AbstractMomentObjective
	sense::MOI.OptimizationSense
	obj::MomExpr{PT}
end

function MomObj(sense::MOI.OptimizationSense, obj::Mom{PT}) where PT<:MT
    return MomObj(sense,convert(MomExpr{PT},obj))
end

function MomObj(sense::MOI.OptimizationSense, pol::PT, meas::Measure) where PT<:MT
    return MomObj(sense,Mom(pol,meas))
end

function MomObj(sense::MOI.OptimizationSense, meas::Measure, pol::PT) where PT<:Union{Number,AbstractPolynomialLike}
	return MomObj(sense,Mom(pol,meas))
end

function Base.show(io::IO,f::MO) where MO<:AbstractMomentObjective
    if typeof(f) == EmptyObjective
        print(io,"Feasability problem:")
    elseif f.sense == MOI.MAX_SENSE
        print(io,"Maximize $(f.obj)")
    else
        print(io,"Minimize $(f.obj)")
    end
end

function measures(f::MO) where MO<:AbstractMomentObjective
	return collect(keys(f.obj.momdict))
end


"""
NLMomObj
"""
#TODO: mutable struct NLMomObj <:AbstractMomentObjective end

"""
MomCon
"""
mutable struct MomCon{PT<:MT} <: JuMP.AbstractConstraint
	func::MomExpr{PT}
    set::MOI.AbstractScalarSet
end

function MomCon(func::AffMomExpr, set::MOI.AbstractScalarSet)
    return MomCon(momexpr(func), MOI.Utilities.shift_constant(set, -constant(func)))
end

function MomCon(func::Mom, set::MOI.AbstractScalarSet)
    return MomCon(MomExpr(func), set)
end

function Base.convert(::Type{MomCon{PT1}},mc::MomCon) where {PT1<:MT}
    return MomCon{PT1}(convert(MomExpr{PT1},mc.func),mc.set)
end

function Base.promote_rule(::Type{MomCon{PT1}},::Type{MomCon{PT2}})  where {PT1<:MT, PT2<:MT}
    return MomCon{promote_type{PT1,PT2}}
end

JuMP.jump_function(con::MomCon) = con.func
JuMP.moi_set(con::MomCon) = con.set
JuMP.shape(con::MomCon) = MomConShape()
JuMP.reshape_vector(expr::MomExpr, ::MomConShape) = expr
JuMP.reshape_set(set::MOI.AbstractScalarSet, ::MomConShape) = set

function Base.show(io::IO, mc::MomCon)
    print(io, JuMP.constraint_string(JuMP.REPLMode, mc))
end

function JuMP.function_string(mode, mc::MomCon)
    return string(mc.func)
end

function measures(mc::MomCon)
    return measures(mc.func.momdict)
end

function measures(mcv::Vector{MomCon})
    measv = Set{Measure}()
	for mc in mcv
        union!(measv, measures(mc))
	end
    return measv
end

#TODO: remove when fixed in MOI
Base.broadcastable(set::MOI.AbstractScalarSet) = Ref(set)
