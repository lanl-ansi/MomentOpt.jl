export MomObj, MomCon, MomCons
# export NLMomObj

abstract type AbstractMomentObjective end
struct EmptyObjective <: AbstractMomentObjective end

abstract type AbstractMomentConstraint end


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
mutable struct MomCon{PT<:MT,T<:Number} <:AbstractMomentConstraint
	lhs::MomExpr{PT}
	mode::Symbol
	rhs::T
end

function MomCon(lhs::Mom{PT}, mode::Symbol, rhs::T) where {T<:Number, PT<:MT}
	return MomCon{PT,T}(convert(MomExpr{PT},lhs),mode,rhs)
end

function MomCon(meas::Measure,pol::PT,mode::Symbol,rhs::T) where {T<:Number,PT<:MT}
	return MomCon(Mom(meas,pol),mode,rhs)
end

function  MomCon(pol::PT,meas::Measure,mode::Symbol,rhs::T) where {T<:Number,PT<:MT}
	return MomCon(Mom(meas,pol),mode,rhs)
end

function Base.convert(::Type{MomCon{PT1,T1}},mc::MomCon) where {T1<:Number, PT1<:MT}
	return MomCon{PT1,T1}(convert(MomExpr{PT1},mc.lhs),mc.mode,convert(T1,mc.rhs))
end


function Base.promote_rule(::Type{MomCon{PT1,T1}},::Type{MomCon{PT2,T2}})  where {T1<:Number,T2<:Number, PT1<:MT, PT2<:MT}
	return MomCon{promote_type{PT1,PT2},promote_type{T1,T2}}
end

function Base.show(io::IO,mc::MC) where MC <: AbstractMomentConstraint
	if mc.mode == :eq
		print(io,"$(mc.lhs) = $(mc.rhs)")
	elseif mc.mode ==:leq
		print(io,"$(mc.lhs) ⩽  $(mc.rhs)")
	elseif mc.mode ==:geq
		print(io,"$(mc.lhs) ⩾ $(mc.rhs)")
	end
end

"""
MomCons
"""

function MomCons(lhs::M,mode::Symbol,rhs::T) where {M<:AbstractMomentExpressionLike, T<:Number}
	return MomCon(lhs,mode,rhs)
end

function MomCons(lhs::Vector{M},mode::Vector{Symbol},rhs::Vector{T}) where {M<:AbstractMomentExpressionLike, T<:Number}
	MT = promote_type([montype(l) for l in lhs]...)
	mc = MomCon{MT,T}[]
	for i in 1:length(rhs)
		push!(mc,MomCon(lhs[i],mode[i],rhs[i]))
	end
	return mc
end

function MomCons(lhs::Vector{M},mode::Symbol,rhs::Vector{T}) where {M<:AbstractMomentExpressionLike, T<:Number}
	MT = promote_type([montype(l) for l in lhs]...)
	mc = MomCon{MT,T}[]
	for i in 1:length(rhs)
		append!(mc,[MomCon(lhs[i],mode,rhs[i])])
	end
	return mc
end

function MomCons(lhs::Vector{M},mode::Symbol,rhs::T) where {M<:AbstractMomentExpressionLike, T<:Number}
 	MT = promote_type([montype(l) for l in lhs]...)
	mc = MomCon{MT,T}[]
	for i in 1:length(rhs)
		append!(mc,[MomCon(lhs[i],mode,rhs)])
	end
	return mc
end


# pretty printing
function Base.show(io::IO,momcons::Vector{M}) where M<:AbstractMomentConstraint
	for i=1:length(momcons)
		println(io, "$(momcons[i])")
	end
end

function measures(mcv::Vector{M}) where M<:AbstractMomentConstraint
	measv = Measure[]
	for mc in mcv
		append!(measv,measures(mc.lhs))
	end
	unique!(measv)
end
function measures(mcv::M) where M<:AbstractMomentConstraint
	measures([mcv])
end
