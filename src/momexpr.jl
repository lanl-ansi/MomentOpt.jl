export Mom, MomExpr, AffMomExpr

abstract type AbstractMomentExpressionLike end
abstract type AbstractMomentExpression <: AbstractMomentExpressionLike end

Base.broadcastable(mom::AbstractMomentExpressionLike) = Ref(mom)

function compatible(meas::Measure,mon::MT)
    if typeof(mon)<:Number
        return true
    elseif isempty(setdiff(variables(mon),variables(meas)))
        return true
    else
        return false
    end
end

"""
    Mom

Type to represent a moment, i.e. the pairing of a measure with a polynomial.
"""
mutable struct Mom{PT<:MT} <: AbstractMomentExpression
    meas::Measure
    mon::PT
    function Mom(meas::Measure,mon::MT)
        if compatible(meas, mon)
            return new{typeof(mon)}(meas,mon)
        else
            @error "$meas does not act on $mon"
        end
    end
end

function Mom(mon)
    @error "The definition of a moment requires a measure and an integrand."
end

function Mom(mon::PT,meas::Measure) where PT<:MT
    return Mom(meas,mon)
end

function Base.dot(mon::PT, meas:Measure) where PT<:MT
    return Mom(meas, mom)
end

function Base.dor(meas::Measure, mon::PT) where PT <: MT
    return Mom(meas, mom)
end

# conversion and promotion
function montype(m::Mom)
    return typeof(m.mon)
end

function Base.copy(m::Mom)
    return Mom(m.meas,copy(m.mon))
end

function Base.convert(::Type{Mom{PT}}, m::Mom) where PT<:MT
    return Mom(m.meas,convert(PT,m.mon))
end

function Base.promote_rule(::Type{Mom{PT1}},::Type{Mom{PT2}}) where {PT1<:MT, PT2<:MT}
    return Mom{promote_type(PT1,PT2)}
end

# pretty printing
function Base.show(io::IO,mom::Mom)
    print(io, "⟨$(mom.meas), $(mom.mon)⟩")
end

"""
    MomExpr

Type to represent a linear combination of Mom.
"""
mutable struct MomExpr{PT<:MT} <: AbstractMomentExpression
    momdict::OrderedDict{Measure,PT}
    function MomExpr(momdict::OrderedDict{Measure,PT}) where PT<:MT
        if all(args->compatible(args...),momdict)
            new{PT}(momdict)
        else
            @error "Measures and moments are not compatible."
        end
    end
end

function MomExpr(momdict::OrderedDict{<:Measure,PT}) where PT<:MT
    return MomExpr(convert(OrderedDict{Measure,PT},momdict))
end

# backwards compatibility of constructor
function MomExpr(poly::PT, mu::Measure) where PT<:MT
    return MomExpr(OrderedDict{Measure,PT}(mu => poly))
end

function MomExpr(mu::Measure,poly::PT) where PT<:MT
    return MomExpr(OrderedDict{Measure,PT}(mu => poly))
end

function MomExpr(mom::Mom)    
    return MomExpr(OrderedDict{Measure,montype(mom)}(mom.meas => mom.mon))
end

# conversion and promotion
function montype(m::ME) where ME <:AbstractMomentExpression
    return promote_type([typeof(m.momdict[meas]) for meas in keys(m.momdict)]...)
end

function Base.convert(::Type{MomExpr{PT}},m::Mom) where PT<:MT
    return MomExpr(OrderedDict{Measure,PT}(m.meas => convert(PT,m.mon)))
end

function Base.convert(::Type{MomExpr{PT}}, m::AbstractMomentExpression) where PT<:MT
    momdict = OrderedDict{Measure,PT}()
    for meas in keys(m.momdict)
        momdict[meas] = convert(PT,m.momdict[meas])
    end
    return MomExpr(momdict)
end

function Base.promote_rule(::Type{MomExpr{PT1}},::Type{Mom{PT2}}) where {PT1<:MT, PT2<:MT}
    return MomExpr{promote_type(PT1,PT2)}
end

# pretty printing
function Base.show(io::IO,me::MomExpr) 
    n = length(me.momdict)
    for m in keys(me.momdict)
        print(io, Mom(m,me.momdict[m]))
        n = n-1
        if n>0
            print(io, " + ")
        end
    end
end

function measures(me::MomExpr)
    return collect(keys(me.momdict))
end

"""
    AffMomExpr

Type for affine moment expressions.
"""
mutable struct AffMomExpr{PT<:MT,T<:Number}  <: AbstractMomentExpressionLike
    exp::MomExpr{PT}
    con::T
end

function momexpr(ae::AffMomExpr)
    return ae.exp
end

function constant(ae::AffMomExpr)
    return ae.con
end

function AffMomExpr(mom::Union{Mom,Number},c::Number)
    return AffMomExpr(MomExpr(mom),c)
end

# conversion and promotion
function Base.promote_rule(::Type{AffMomExpr{PT1,T1}},::Type{AffMomExpr{PT2,T2}}) where {PT1<:MT, PT2<:MT, T1<:Number, T2<:Number}
    return AffMomExpr{promote_type(PT1,PT2), promote_type(T1,T2)}
end

function Base.convert(::Type{AffMomExpr{PT,T}}, ae::AffMomExpr) where {PT<:MT, T<:Number}
    return AffMomExpr(convert(MomExpr{PT}, momexpr(ae)), convert(T,constant(ae)))
end

function Base.promote_rule(::Type{AffMomExpr{PT1,T1}},::Type{MomExpr{PT2}}) where {PT1<:MT, PT2<:MT, T1<:Number}
    return AffMomExpr{promote_type(PT1,PT2), promote_type(T1,Int)}
end

function Base.convert(::Type{AffMomExpr{PT,T}}, me::MomExpr) where {PT<:MT, T<:Number}
    return AffMomExpr(convert(MomExpr{PT}, me), zero(T))
end

function Base.promote_rule(::Type{AffMomExpr{PT1,T1}},::Type{Mom{PT2}}) where {PT1<:MT, PT2<:MT, T1<:Number}
    return AffMomExpr{promote_type(PT1,PT2), promote_type(T1,Int)}
end

function Base.convert(::Type{AffMomExpr{PT,T}}, mom::Mom) where {PT<:MT, T<:Number}
    return AffMomExpr(convert(MomExpr{PT}, MomExpr(mom)), zero(T))
end

# pretty printing
function Base.show(ae::AffMomExpr)
    if constant(ae)>0
        print(io,"$(momexpr(ae)) + $(constant(ae))")
    elseif constant(ae)<0
        print(io,"$(momexpr(ae)) - $(-(constant(ae)))")
    else
        print(io, momexpr(ae))
    end
end

"""
Linear operations
"""

# Mom
function Base.:*(a::Number, m::Mom)
    return Mom(m.meas,a*copy(m.mon))
end

function Base.:-(mom::Mom)
    return Mom(mom.meas,-mom.mon)
end

# MomExpr
function Base.:*(a::Number, me::MomExpr{PT}) where PT
    PU = typeof(one(typeof(a)) * one(PT))
    return MomExpr(OrderedDict{Measure, PU}(meas => a*poly for (meas,poly) in me.momdict))
end

function Base.:-(me::MomExpr{PT}) where PT
    PU = typeof(-one(PT))
    return MomExpr(OrderedDict{Measure, PU}(meas => -poly for (meas, poly) in me.momdict))
end

# AffMomExpr
function Base.:*(a::Number, ae::AffMomExpr)
    return AffMomExpr(a*momexpr(ae), a*constant(ae))
end

function Base.:-(ae::AffMomExpr)
    return AffMomExpr(-momexpr(ae), -constant(ae))
end

# AbstractMomentExpressionLike
function Base.:*(amel::AbstractMomentExpressionLike, a::Number)
    return a*amel
end

function Base.:/(amel::AbstractMomentExpressionLike, a::Number)
    return (1/a)*amel
end


# summation
function add_mom_type(mev::Vector{<:Union{Mom{T}, MomExpr{T}}}) where T<:Number
    return T
end

function add_mom_type(mev::Vector{<:Union{Mom{T}, MomExpr{T}}}) where T<:AbstractPolynomialLike
    return polynomialtype(T)
end

function Base.sum(mev::Vector{<:Mom})
    T = add_mom_type(mev)
    momdict = OrderedDict{Measure,T}()
    for me in mev
        if haskey(momdict,me.meas)
            momdict[me.meas] = momdict[me.meas] + me.mon
        else
            momdict[me.meas] = me.mon
        end
    end
    return MomExpr(momdict)
end

function Base.sum(mev::Vector{<:MomExpr})
    T = add_mom_type(mev)
    momdict = OrderedDict{Measure,T}()
    for me in mev
        merge!(+, momdict, me.momdict)
    end
    return MomExpr(momdict)
end

function Base.sum(aev::Vector{<:AffMomExpr})
    return AffMomExpr(sum(momexpr.(aev)), sum(constant.(aev)))
end

# plus
function Base.:+(mom1::AbstractMomentExpressionLike,mom2::AbstractMomentExpressionLike)
    return sum([mom1,mom2])
end

function Base.:+(me::Mom, c::Number)
    return AffMomExpr(MomExpr(me),c)
end

function Base.:+(me::MomExpr, c::Number)
    return AffMomExpr(me,c)
end

function Base.:+(me::AffMomExpr, c::Number)
    return AffMomExpr(momexpr(me),c+constant(me))
end

function Base.:+(c::Number,ae::AbstractMomentExpressionLike)
    return ae+c
end

# minus
function Base.:-(ae1::AbstractMomentExpressionLike,ae2::Union{Number,AbstractMomentExpressionLike})
    return ae1+(-ae2)
end

function Base.:-(ae1::Union{Number,AbstractMomentExpressionLike},ae2::AbstractMomentExpressionLike)
    return ae1+(-ae2)
end

function Base.:-(ae1::AbstractMomentExpressionLike,ae2::AbstractMomentExpressionLike)
    return ae1+(-ae2)
end
