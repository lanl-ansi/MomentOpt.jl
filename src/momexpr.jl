export CMom, Mom, CMomExpr, MomExpr

abstract type AbstractMomentExpressionLike end
abstract type AbstractMomentExpression <: AbstractMomentExpressionLike end
abstract type AbstractMoment  <: AbstractMomentExpressionLike end

"""
CMom 
"""
mutable struct CMom{T <: Number} <: AbstractMoment
	meas::Measure
       	mon::T 
end

function CMom(mon::T,meas::Measure) where T<: Number
	return Cmom(meas,mon)
end
# constructor for vectors
function CMom(meas::Measure,monv::Vector{T}) where T <: Number
	monv = Vector{CMom{T}}()
	for mon in monv
		push!(momv, CMom(meas,mon))
	end
	return momv
end

function CMom(monv::Vector{T},meas::Measure) where T <: Number
	return Mom(meas,monv)
end


# conversion and promotion
function Base.convert(::Type{CMom{T1}},cmon::CMom) where {T1<:Number}
	return CMom(cmon.meas,convert(T1,cmon.mon))
end

function Base.promote_rule(::Type{CMom{T1}},::Type{CMom{T2}}) where {T1<:Number, T2<:Number}
	return CMom{promote_type(T1,T2)}
end

"""
Mom
"""
function compatible(meas::Measure,mon::PT) where PT<:Union{Number,AbstractPolynomialLike}
	if typeof(mon)<:Number
		return true
	elseif isempty(setdiff(MP.variables(mon),variables(meas)))
		return true
	else
		return false
	end
end
	
mutable struct Mom{PT<:AbstractPolynomialLike} <: AbstractMoment
	meas::Measure
	mon::PT
	function Mom(meas::Measure,mon::PT) where PT<:AbstractPolynomialLike
		if compatible(meas, mon)
			return new{PT}(meas,mon)
		else
			@error "$meas does not act on $(mon)"
		end
	end
end

function Mom(mon::PT,meas::Measure) where PT<:AbstractPolynomialLike
	return Mom(meas,mon)
end

# backwards compatibility of constructor
function Mom(meas::Measure, c::T) where T<:Number
	return CMom(meas,c)
end

function Mom(c::T,meas::Measure) where T<:Number
	return CMom(meas,c)
end

# constructor for vectors
function Mom(meas::Measure,monv::Vector{PT}) where PT <: Union{AbstractPolynomialLike,Number}
	momv = Vector{Mom{PT}}()
	for mon in monv
		push!(momv, Mom(meas,mon))
	end
	return momv
end

function Mom(monv::Vector{PT},meas::Measure) where PT <: Union{AbstractPolynomialLike,Number}
	return Mom(meas,monv)
end

# conversion and promotion
function montype(m::MOM) where MOM <: AbstractMoment
	return typeof(m.mon)
end
 
function Base.convert(::Type{Mom{PT}}, m::M) where {PT<:AbstractPolynomialLike,M<:AbstractMoment}
	return Mom(m.meas,convert(PT,m.mon))
end

function Base.promote_rule(::Type{Mom{PT}},::Type{CMom{T}}) where {T<:Number, PT<: AbstractPolynomialLike}
	return Mom{promote_type(PT,T)}
end

function Base.promote_rule(::Type{Mom{PT1}},::Type{Mom{PT2}}) where {PT1<:AbstractPolynomialLike, PT2<: AbstractPolynomialLike}
	return Mom{promote_type(PT1,PT2)}
end


# abstract type stable arithmetic
function Base.:*(a::T, m::MOM) where {T<:Number,MOM<:AbstractMoment}
	return Mom(m.meas,a*copy(m.mon))
end

function Base.:*(m::MOM,a::T)  where {T<:Number,MOM<:AbstractMoment}
	return Mom(m.meas,a*copy(m.mon))
end

function Base.:/(m::MOM, a::T) where {T<:Number,MOM<:AbstractMoment}
	return Mom(m.meas,copy(m.mon)/a)
end

function NegMom(mom::MOM) where MOM <: AbstractMoment
	return Mom(mom.meas,-mom.mon)
end

# pretty printing
function Base.show(io::IO,mom::MOM) where MOM<:AbstractMoment
	print(io, "< $(mom.meas), $(mom.mon) >")
end


"""
CMomExpr 
"""
mutable struct CMomExpr{T<:Number} <: AbstractMomentExpression
	momdict::Dict{Measure,T} 
end

function CMomExpr(poly::T, mu::Measure) where T<:Number
	return CMomExpr(Dict{Measure,T}(mu=>poly))
end

function CMomExpr(mu::Measure,poly::T) where T<:Number
	return CMomExpr(Dict{Measure,T}(mu=>poly))
end

function CMomExpr(cmom::CMom{T}) where T<:Number
	return CMomExpr(cmom.meas,cmom.mon)
end

# constructors for vectors
function CMomExpr(meas::Measure,cv::Vector{T}) where T<:Number
	mev = CMomExpr{T}[]
	for c in cv
		push!(mev, CMomExpr(meas,c))
	end
	return mev 
end

function CMomExpr(cv::Vector{T},meas::Measure) where T<: Number
	return CMomExpr(meas,cv)
end

function CMomExpr(mom::Vector{CMom{T}}) where T<:Number
	mev = CMomExpr{T}[]        
	for m in mom
		push!(mev, CMomExpr(m))
	end
	return mev 
end


# conversion and promotion
function Base.convert(::Type{CMomExpr{T}},cm::CMom) where T<:Number
	return CMomExpr(Dict{Measure,T}(cm.meas=>convert(T,cm.mon)))
end

function Base.convert(::Type{CMomExpr{T1}}, cme::CMomExpr) where T1<: Number
	ncme = Dict{Measure,T1}()
	for meas in keys(cme.momdict)
		ncme[meas]= convert(T1,cme.momdict[meas])
	end
	return CMomExpr(ncme)
end

function Base.promote_rule(::Type{CMomExpr{T1}},::Type{CMom{T2}}) where {T1<: Number, T2<: Number}
	return CMomExpr{promote_type(T1,T2)}
end

function Base.promote_rule(::Type{CMomExpr{T1}},::Type{CMomExpr{T2}}) where {T1<: Number, T2<: Number}
	return CMomExpr{promote_type(T1,T2)}
end



"""
MomExpr 
"""
mutable struct MomExpr{PT<:AbstractPolynomialLike}  <: AbstractMomentExpression
	momdict::Dict{Measure,PT}
	function MomExpr(momdict::Dict{Measure,PT}) where PT<:AbstractPolynomialLike
		bool = true
		for meas in keys(momdict)
		bool = bool && compatible(meas,momdict[meas])
		end
		if bool
		new{PT}(momdict)
		else
		@error "measures and moments are not compatible"
		return nothing
		end
	end
end

# backwards compatibility of constructor
function MomExpr(momdict::Dict{Measure,T}) where T<: Number
	return CMomExpr(momdict)
end

function MomExpr(poly::PT, mu::Measure) where PT<:Union{AbstractPolynomialLike,Number}
	return MomExpr(Dict{Measure,PT}(mu=>poly))
end

function MomExpr(mu::Measure,poly::PT) where PT<:Union{AbstractPolynomialLike,Number}
	return MomExpr(Dict{Measure,PT}(mu=>poly))
end

function MomExpr(mom::MOM) where MOM<:AbstractMoment
	return MomExpr(Dict{Measure,typeof(mom.mon)}(mom.meas=>mom.mon))
end


function montype(m::MOM) where MOM <: AbstractMoment
	return typeof(m.mon)
end

# constructors for vectors
function MomExpr(mom::Vector{MOM}) where MOM<:AbstractMoment
	MT = promote_type([montype(m) for m in mom]...)
	mev = MomExpr{MT}[]
	for m in mom
		push!(mev, MomExpr(m))
	end
	return mev 
end
function MomExpr(meas::Measure,monv::Vector{PT}) where PT <: Union{AbstractPolynomialLike,Number}
	MT = promote_type([montype(m) for m in mom]...)
	mev = MomExpr{MT}()
	for m in mom
		push!(momv, Mom(meas,mon))
	end
	return momv
end

function MomExpr(monv::Vector{PT},meas::Measure) where PT <: Union{AbstractPolynomialLike,Number}
	return Mom(meas,monv)
end


# conversion and promotion
function montype(m::ME) where ME <:AbstractMomentExpression
	return promote_type([typeof(m.momdict[meas]) for meas in keys(m.momdict)]...)
end

function Base.convert(::Type{MomExpr{PT}},m::MOM) where {PT<:AbstractPolynomialLike, MOM<:AbstractMoment}
	return MomExpr(Dict{Measure,PT}(m.meas=> convert(PT,m.mon)))
end

function Base.convert(::Type{MomExpr{PT}}, m::ME) where {PT<:AbstractPolynomialLike,ME<: AbstractMomentExpression}
	momdict = Dict{Measure,PT}()
	for meas in keys(m.momdict)
		momdict[meas] = convert(PT,m.momdict[meas])
	end
	return MomExpr(momdict)
end

function Base.promote_rule(::Type{CMomExpr{T1}},::Type{CMom{T2}}) where {T1<:Number,T2<:Number}
	return CMomExpr{promote_type(T1,T2)}
end

function Base.promote_rule(::Type{MomExpr{PT}},::Type{CMom{T}}) where {T<:Number, PT<: AbstractPolynomialLike}
	return MomExpr{promote_type(PT,T)}
end
function Base.promote_rule(::Type{MomExpr{PT1}},::Type{Mom{PT2}}) where {PT1<:AbstractPolynomialLike, PT2<:AbstractPolynomialLike}
	return MomExpr{promote_type(PT1,PT2)}
end

function Base.promote_rule(::Type{MomExpr{PT}},::Type{CMomExpr{T}}) where {T<:Number, PT<: AbstractPolynomialLike}
	return MomExpr{promote_type(PT,T)}
end

function Base.promote_rule(::Type{MomExpr{PT1}},::Type{MomExpr{PT2}}) where {PT1<:AbstractPolynomialLike, PT2<:AbstractPolynomialLike}
	return MomExpr{promote_type(PT1,PT2)}
end



# abstract type stable arithmetic
function Base.:*(a::T, me::ME) where {T<:Number, ME<:AbstractMomentExpression}
	nme = MomExpr(copy(me.momdict))
	for (meas,poly) in me.momdict
		nme.momdict[meas] = a*poly
	end
	return nme
end

function Base.:*(me::ME,a::T) where  {T<:Number, ME<:AbstractMomentExpression}
	return a*me 
end

function Base.:/(me::ME,a::T) where {T<:Number, ME<:AbstractMomentExpression} 
	return (1/a)*me
end

function NegMomExpr(mom::ME) where ME <: AbstractMomentExpression
	mmom = MomExpr(copy(mom.momdict))
	for meas in keys(mom.momdict)
		mmom.momdict[meas] = -mom.momdict[meas]
	end
	return mmom
end

function NegMomExpr(mom::MOM) where MOM <:AbstractMoment
	return NegMom(mom)
end


# pretty printing
function Base.show(io::IO,me::ME) where ME<:AbstractMomentExpression
	n = length(me.momdict)
	for m in keys(me.momdict)
		print(io, "< $m, $(me.momdict[m]) >")
		n = n-1
		if n>0
		print(io, " + ")
		end
	end
end

function measures(me::ME) where ME<:AbstractMomentExpression
	return collect(keys(me.momdict))
end


"""
Type unstabil operations
"""

function AddMom(mev::Vector{MOM}) where MOM<:AbstractMoment
	T = typeof(sum([me.mon for me in mev]))
	momdict = Dict{Measure,T}()
	for me in mev
		if haskey(momdict,me.meas)
		momdict[me.meas] = momdict[me.meas] + me.mon
		else
		momdict[me.meas] = me.mon
		end
	end
	if momdict.count ==1
		for meas in keys(momdict)
		return Mom(meas, momdict[meas])
		end
	else
		return MomExpr(momdict)
	end
end

function AddMomExpr(mev::Vector{ME}) where ME <: AbstractMomentExpression
	T = promote_type([typeof(me.momdict) for me in mev]...)
	momdict = T()
	for me in mev
		for meas in keys(me.momdict)
			if haskey(momdict,meas)
			momdict[meas] = momdict[meas] + me.momdict[meas]
			else
			momdict[meas] = me.momdict[meas]
			end
		end
	end
	return MomExpr(momdict)
end

function AddMomExpr(mev::Vector{MOM}) where MOM<: AbstractMoment
	return AddMom(mev)
end

function Base.:+(mom1::MEL1,mom2::MEL2) where {MEL1<:AbstractMomentExpressionLike, MEL2<:AbstractMomentExpressionLike}
	return AddMomExpr([mom1,mom2])
end

function Base.:-(me1::MEL1,me2::MEL2) where {MEL1<:AbstractMomentExpressionLike,MEL2<:AbstractMomentExpressionLike}
	return me1+NegMomExpr(me2)
end




# Vector operations


function Base.:+(mv1::Vector{MEL1},mv2::Vector{MEL2}) where {MEL1<:AbstractMomentExpressionLike,MEL2<:AbstractMomentExpressionLike}
	if length(mv1)==length(mv2)
	T = promote_type([montype(m) for m in mv1]...,[montype(m) for m in mv2]...)
	MT = Base.promote_op(+,T,T)
	mv = MomExpr{MT}[]
	for i=1:length(mv1)
		push!(mv, mv1[i]+mv2[i])
	end
	return mv
	else
		println("Vectors need to have same length")
		return nothing
	end
end

function Base.:-(mv1::Vector{MEL1},mv2::Vector{MEL2}) where {MEL1<:AbstractMomentExpressionLike,MEL2<:AbstractMomentExpressionLike}
	if length(mv1)==length(mv2)
	T = promote_type([montype(m) for m in mv1]...,[montype(m) for m in mv2]...)
	MT = Base.promote_op(+,T,T)
	mv = MomExpr{MT}[]
	for i=1:length(mv1)
		push!(mv, mv1[i]-mv2[i])
	end
	return mv
	else
		println("Vectors need to have same length")
		return nothing
	end
end

function Base.:*(mv::Vector{MEL},a::T) where {T<:Number, MEL<:AbstractMomentExpressionLike}
	NT = promote_type(MEL,CMom{T})
	nmv= NT[]
	for i=1:length(mv)
		push!(nmv, a*mv[i])
	end
	return mv
end

function Base.:*(a::T,mv::Vector{MEL}) where {T<:Number, MEL<:AbstractMomentExpressionLike}
	return mv*a
end
function Base.:/(mv::Vector{MEL},a::T) where {T<:Number, MEL<:AbstractMomentExpressionLike}
	return mv*(1/a)
end


