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


# abstract type stable arithmetic
function Base.:*(a::Number, m::Mom)
    return Mom(m.meas,a*copy(m.mon))
end

function Base.:*(m::Mom,a::Number) 
    return Mom(m.meas,a*copy(m.mon))
end

function Base.:/(m::Mom, a::Number) 
    return Mom(m.meas,copy(m.mon)/a)
end

function Base.:-(mom::Mom)
    return Mom(mom.meas,-mom.mon)
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
    momdict::Dict{Measure,PT}
    function MomExpr(momdict::Dict{Measure,PT}) where PT<:MT
        if all(args->compatible(args...),momdict)
            new{PT}(momdict)
        else
            @error "Measures and moments are not compatible."
        end
    end
end

function MomExpr(momdict::Dict{<:Measure,PT}) where PT<:MT
    return MomExpr(convert(Dict{Measure,PT},momdict))
end

# backwards compatibility of constructor
function MomExpr(poly::PT, mu::Measure) where PT<:MT
    return MomExpr(Dict{Measure,PT}(mu => poly))
end

function MomExpr(mu::Measure,poly::PT) where PT<:MT
    return MomExpr(Dict{Measure,PT}(mu => poly))
end

function MomExpr(mom::Mom)    
    return MomExpr(Dict{Measure,montype(mom)}(mom.meas => mom.mon))
end

# conversion and promotion
function montype(m::ME) where ME <:AbstractMomentExpression
    return promote_type([typeof(m.momdict[meas]) for meas in keys(m.momdict)]...)
end

function Base.convert(::Type{MomExpr{PT}},m::Mom) where PT<:MT
    return MomExpr(Dict{Measure,PT}(m.meas => convert(PT,m.mon)))
end

function Base.convert(::Type{MomExpr{PT}}, m::AbstractMomentExpression) where PT<:MT
    momdict = Dict{Measure,PT}()
    for meas in keys(m.momdict)
        momdict[meas] = convert(PT,m.momdict[meas])
    end
    return MomExpr(momdict)
end

function Base.promote_rule(::Type{MomExpr{PT1}},::Type{Mom{PT2}}) where {PT1<:MT, PT2<:MT}
    return MomExpr{promote_type(PT1,PT2)}
end

# abstract type stable arithmetic
function Base.:*(a::T, me::MomExpr) where {T<:Number}
    return MomExpr(Dict(meas => a*poly for (meas,poly) in me.momdict))
end

function Base.:*(me::MomExpr,a::T) where  {T<:Number}
    return a*me
end

function Base.:/(me::MomExpr,a::T) where {T<:Number}
    return (1/a)*me
end

function Base.:-(mom::MomExpr)
    return MomExpr(Dict(meas => -poly for (meas, poly) in mom.momdict))
end

function add_mom_type(mev::Vector{<:Union{Mom{T},MomExpr{T}}}) where T<:Number
    return T
end

function add_mom_type(mev::Vector{<:Union{Mom{T},MomExpr{T}}}) where T<:AbstractPolynomialLike
    return polynomialtype(T)
end

function Base.sum(mev::Vector{<:Mom})
    T = add_mom_type(mev)
    momdict = Dict{Measure,T}()
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
    momdict = Dict{Measure,T}()
    for me in mev
        merge!(+, momdict, me.momdict)
    end
    return MomExpr(momdict)
end

function Base.:+(mom1::AbstractMomentExpression,mom2::AbstractMomentExpression)
    return sum([mom1,mom2])
end

function Base.:-(mom1::AbstractMomentExpression,mom2::AbstractMomentExpression)
   return mom1+(-mom2)
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
    momexp::MomExpr{PT}
    cons::T
end

function constant(ame::AffMomExpr)
    return ame.cons
end

function momexpr(ame::AffMomExpr)
    return ame.momexp
end
