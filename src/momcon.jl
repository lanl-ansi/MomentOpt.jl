export MomObj, CMomObj, MomCon, CMomCon, MomCons
export EmptyObjective
# export NLMomObj

abstract type AbstractMomentObjective end
struct EmptyObjective <: AbstractMomentObjective end

abstract type AbstractMomentConstraint end


"""
(C)MomObj
"""
mutable struct MomObj{PT<:AbstractPolynomialLike} <:AbstractMomentObjective
	sense::Symbol
	obj::MomExpr{PT}
end

function MomObj(sense::Symbol, obj::Mom{PT}) where PT<:AbstractPolynomialLike
	return MomObj(sense,convert(MomExpr{PT},obj))
end

function MomObj(sense::Symbol, obj::Mom{PT}) where PT<:AbstractPolynomialLike
	return MomObj(sense,convert(MomExpr{PT},obj))
end

function MomObj(sense::Symbol, pol::PT, meas::Measure) where PT<:Union{Number,AbstractPolynomialLike}
	return MomObj(sense,Mom(pol,meas))
end
function MomObj(sense::Symbol, meas::Measure,pol::PT) where PT<:Union{Number,AbstractPolynomialLike}
	return MomObj(sense,Mom(pol,meas))
end

mutable struct CMomObj{T<:Number} <:AbstractMomentObjective
	sense::Symbol
	obj::CMomExpr{T}
end

function MomObj(sense::Symbol, obj::CMom{T}) where T<:Number
	return CMomObj(sense,convert(CMomExpr{T},obj))
end

function CMomObj(sense::Symbol, obj::CMom{PT}) where PT<:AbstractPolynomialLike
	return MomObj(sense,convert(CMomExpr{PT},obj))
end

function Base.show(io::IO,f::MO) where MO<:AbstractMomentObjective
	if typeof(f)==EmptyObjective
	print(io,"Feasability problem:")
	else
	print(io,"$(f.sense) $(f.obj)")
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
(C)MomCon
"""
mutable struct CMomCon{T1<:Number,T2<:Number} <:AbstractMomentConstraint
	lhs::CMomExpr{T1}
	mode::Symbol
	rhs::T2
end

function CMomCon(cmom::CMom,mode::Symbol,rhs::T) where T<: Number
	return CMomCon(CMomExpr(cmom),mode,rhs)
end

function Base.convert(::Type{CMomCon{T1,T2}},mc::CMomCon) where {T1<:Number, T2<:Number}
	return CMomCon{T1,T2}(convert(CMomExpr{T1},mc.lhs),mc.mode,convert(T2,mc.rhs))
end

function Base.promote_rule(::Type{CMomCon{T1,T2}},::Type{CMomCon{T3,T4}}) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}
	return CMomCon{promote_type{T1,T3},promote_type{T2,T4}}
end

mutable struct MomCon{PT<:AbstractPolynomialLike,T<:Number} <:AbstractMomentConstraint
	lhs::MomExpr{PT}
	mode::Symbol
	rhs::T
end

function MomCon(lhs::Mom{PT}, mode::Symbol, rhs::T) where {T<:Number, PT<:AbstractPolynomialLike}
	return MomCon{PT,T}(convert(MomExpr{PT},lhs),mode,rhs)
end

function MomCon(lhs::CMom{T1}, mode::Symbol, rhs::T2) where {T1<:Number, T2<:Number}
	return CMomCon{T1,T2}(convert(CMomExpr{T1},lhs),mode,rhs)
end

function MomCon(meas::Measure,pol::PT,mode::Symbol,rhs::T) where {T<:Number,PT<:Union{AbstractPolynomialLike, Number}}
	return MomCon(Mom(meas,pol),mode,rhs)
end

function  MomCon(pol::PT,meas::Measure,mode::Symbol,rhs::T) where {T<:Number,PT<:Union{AbstractPolynomialLike, Number}}
	return MomCon(Mom(meas,pol),mode,rhs)
end

function Base.convert(::Type{MomCon{PT1,T1}},mc::MomCon) where {T1<:Number, PT1<:AbstractPolynomialLike}
	return MomCon{PT1,T1}(convert(MomExpr{PT1},mc.lhs),mc.mode,convert(T1,mc.rhs))
end

function Base.convert(::Type{MomCon{PT1,T1}},mc::CMomCon) where {T1<:Number, PT1<:AbstractPolynomialLike}
	return MomCon{PT1,T1}(convert(MomExpr{PT1},mc.lhs),mc.mode,convert(T1,mc.rhs))
end

function Base.promote_rule(::Type{MomCon{PT1,T1}},::Type{MomCon{PT2,T2}})  where {T1<:Number,T2<:Number, PT1<:AbstractPolynomialLike, PT2<:AbstractPolynomialLike}
	return MomCon{promote_type{PT1,PT2},promote_type{T1,T2}}
end

function Base.promote_rule(::Type{MomCon{PT,T1}},::Type{CMomCon{T2,T3}}) where {PT<:AbstractPolynomialLike,T1<:Number,T2<:Number, T3<:Number}
	return MomCon{promote_type{PT,T2},promote_type{T1,T3}}
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
