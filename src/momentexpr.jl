export Mom
export MomentExpr, AffineMomentExpr

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
    return isempty(measures(meas)) ? true : compatible(first(measures(meas)), mon)
end

"""
    MomentExpr(integrs::Vector{PT}, mexprs::Vector{MeasureExpr}) 

Represents a moment expression.
"""
mutable struct MomentExpr{T <: Number, PT <: MT} <: AbstractMomentExpressionLike
    integr::Vector{PT}
    mexprs::Vector{MeasureExpr{T}}
    function MomentExpr(integrs::Vector{PT}, mexprs::Vector{MeasuresExpr{T}}) where {PT <: MT, T <: Number}
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
            print(io," + ⟨" p, m, "⟩")
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

function MomentExpr(integr::PT, m::AbstractMeasureLIke) where {PT}
    return MomentExpr(integr, MomentExpr([1], [m]))
end

# Linear operations
function Base.:*(a::Number, m::MomentExpr)
    return MeasureExpr(a.*integrands(m), measures(m))
end

function Base.:*(m::MomentExpr, a::Number)
    return a*m
end

function Base.:/(m::MomentExpr, a::Number)
    return 1/a*m
end

function Base.:-( m::MomentExpr)
    return MomentExpr(-integrands(m), measures(m))
end

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

function Base.:+(m1::MomentExpr, m2::MomentExpr)
    return sum([m1, m2])
end

function Base.:-(m1::MomentExpr, m2::MomentExpr)
    return sum([m1, -m2])
end

"""
    sort_by_measure(m::MomentExpr)

Sorts the moments of a moment expression by measure.
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



"""
    sort_by_monomial(m::MomentExpr)

Sorts the moments of a moment expression by monomials.
"""








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

