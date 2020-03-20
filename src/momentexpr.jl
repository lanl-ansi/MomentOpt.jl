export Mom
export sort_by_measure

abstract type AbstractMomentExpressionLike end
Base.broadcastable(mom::AbstractMomentExpressionLike) = Ref(mom)

function compatible(meas::AbstractMeasureLike, mon::MT)
    if typeof(mon) <: Number
        return true
    elseif isempty(setdiff(variables(mon), poly_variables(meas)))
        return true
    else
        return false
    end
end

function compatible(meas::MeasureExpr, mon::MT)
    return isempty(measures(meas)) ? false : compatible(first(measures(meas)), mon)
end

"""
    MomentExpr(integrs::Vector{PT}, mexprs::Vector{MeasureExpr}) 

Represents a moment expression.
"""
mutable struct MomentExpr{T <: Number, PT <: MT} <: AbstractMomentExpressionLike
    integr::Vector{PT}
    mexprs::Vector{MeasureExpr{T}}
    function MomentExpr(integrs::Vector{PT}, mexprs::Vector{MeasureExpr{T}}) where {PT <: MT, T <: Number}
        @assert length(integrs) == length(mexprs) "Inputs need to have same length."
        @assert all(compatible(mexprs[i], integrs[i]) for i = 1:length(integrs)) "Invalid integration, variables do not coincide."
        return new{T, PT}(integrs, mexprs)
    end
end

integrands(me::MomentExpr) = me.integrs
measures(me::MomentExpr) = me.mexprs

function Base.show(io::IO, me::MomentExpr)
    if isempty(integrands(me))
        print(io, 0)
    else
        m = measures(me)[1]
        p = integrands(me)[1]
        print(io,"⟨", p, m, "⟩")
        for (c, m) in zip(integrands(me)[2:end], measures(me)[2:end])
            print(io," + ⟨", p, m, "⟩")
        end
    end
end

# compatibility

function Base.promote_rule(::Type{MomentExpr{T1, PT1}}, ::Type{MomentExpr{T2, PT2}}) where {T1, T2, PT1, PT2}
    return MomenExpr{promote_type(T1, T2), promote_type(PT1, PT2)}
end

function Base.convert(::Type{MomentExpr{T, PT}}, m::MomentExpr) where {T, PT}
    return MomentExpr(convert.(PT, integrands(m)), convert.(MeasureExpr{T}, measures(m)))    
end

function MomentExpr(integr::PT, mexpr::MeasureExpr{T}) where {T, PT}
    return MomentExpr([integr], [mexpr])
end

function MomentExpr(integr::PT, m::AbstractMeasureLike) where {PT}
    return MomentExpr(integr, MomentExpr([1], [m]))
end

"""
    Mom(arg1, arg2)

Construct a MomentExpression.
"""
function Mom(arg1, arg2)
    if typeof(arg1) isa Union{MT, Vector{MT}}
        return MomentExpr(arg1, arg2)
    else
        return MomentExpr(arg2, arg1)
    end
end


"""
    AffineMomentExpr

Type for affine moment expressions.
"""
mutable struct AffineMomentExpr{T <: Number, PT <: MT, S <: Number}  <: AbstractMomentExpressionLike
    exp::MomentExpr{T, PT}
    con::S
end

function momexpr(ae::AffineMomentExpr)
    return ae.exp
end

function constant(ae::AffineMomentExpr)
    return ae.con
end

# conversion and promotion
    function Base.promote_rule(::Type{AffineMomentExpr{T1, PT1, S1}}, ::Type{AffMomExpr{T2, PT2, S2}}) where {T1<:Number, T2<:Number, PT1<:MT, PT2<:MT, S1<:Number, S2<:Number}
    return AffineMomentExpr{promote_type(T1,T2), promote_type(PT1,PT2), promote_type(S1,S2)}
end

function Base.convert(::Type{AffineMomentExpr{T, PT, S}}, ae::AffineMomentExpr) where {T <: Number, PT <: MT, S <: Number}
    return AffineMomentExpr(convert(MomentExpr{T, PT}, momexpr(ae)), convert(s,constant(ae)))
end

function Base.promote_rule(::Type{AffineMomentExpr{T1, PT1, S}},::Type{MomExpr{T2, PT2}}) where {T1 <: Number, T2 <: Number, PT1 <: MT, PT2 <: MT, S <: Number}
    return AffineMomentExpr{promote_type(T1,T2),promote_type(PT1, PT2), promote_type(S,Int)}
end

function Base.convert(::Type{AffineMomentExpr{T, PT, S}}, me::MomExpr) where {T <: Number, PT <: MT, S <: Number}
    return AffineMomentExpr(convert(MomentExpr{T, PT}, me), zero(S))
end

# pretty printing
function Base.show(io::IO, ae::AffineMomentExpr)
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

# MomentExpr

function Base.:*(a::Number, m::MomentExpr)
    return MeasureExpr(a.*integrands(m), measures(m))
end

function Base.:-( m::MomentExpr)
    return MomentExpr(-integrands(m), measures(m))
end

# AffineMomentExpr
function Base.:*(a::Number, ae::AffineMomentExpr)
    return AffineMomentExpr(a*momexpr(ae), a*constant(ae))
end

function Base.:-(ae::AffineMomentExpr)
    return AffineMomentExpr(-momexpr(ae), -constant(ae))
end

# AbstractMomentExpressionLike
function Base.:*(amel::AbstractMomentExpressionLike, a::Number)
    return a*amel
end

function Base.:/(amel::AbstractMomentExpressionLike, a::Number)
    return (1/a)*amel
end


# summation
function Base.sum(mev::Vector{<:MomentExpr})
    integrs = copy(integrands(mev[1]))
    meass = copy(measures(mev[1]))
    for me in mev[2:end]
        for (c, m) in zip(integrands(me), measures(me))
            index = findfirst(==(m), meass)
            if isempty(index)
                push!(integrs, c)
                push!(meass, m)
            else
                integrs[index] += c
            end
        end
    end
    index = findall( x -> !(x==0.0), integrs)
    return MeasureExpr(integrs[index], meass[index])
end

function Base.sum(aev::Vector{<:AffineMomentExpr})
    return AffineMomentExpr(sum(momexpr.(aev)), sum(constant.(aev)))
end

# plus
function Base.:+(mom1::AbstractMomentExpressionLike, mom2::AbstractMomentExpressionLike)
    return sum([mom1, mom2])
end

function Base.:+(me::MomentExpr, c::Number)
    return AffineMomentExpr(me,c)
end

function Base.:+(me::AffineMomentExpr, c::Number)
    return AffineMomentExpr(momexpr(me), c + constant(me))
end

function Base.:+(c::Number, ae::AbstractMomentExpressionLike)
    return ae + c
end

# minus
function Base.:-(ae1::AbstractMomentExpressionLike, ae2::Union{Number, AbstractMomentExpressionLike})
    return ae1 + (-ae2)
end

function Base.:-(ae1::Union{Number, AbstractMomentExpressionLike}, ae2::AbstractMomentExpressionLike)
    return ae1 + (-ae2)
end

function Base.:-(ae1::AbstractMomentExpressionLike, ae2::AbstractMomentExpressionLike)
    return ae1 + (-ae2)
end

"""
    sort_by_measure(m::AbstractMomentExprLike)

Sorts the moments of an (affine) moment expression by measure.
"""
function sort_by_measure(m::MomentExpr{T, PT}) where {T, PT}
    mev = MomentExpr{T, PT}[]
    for (p, mv) in zip(integrands(m), measures(m))
        for (c, meas) in mv
            push!(mev, MomentExpr(v*p, meas))
        end
    end
    return sum(mev)
end

function sort_by_measure(m::AffineMomentExpr)
    return AffineMomentExpr(sort_by_measure(momexpr(m)), constant(m))
end
