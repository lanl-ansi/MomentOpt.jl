export MeasureExpr, AffineMeasureExpr

abstract type AbstractMeasureExpressionLike end
Base.broadcastable(me::AbstractMeasureExpressionLike) = Ref(me)
mutable struct MeasureExpr{T <: Number} <: AbstractMeasureExpressionLike
    coefs::Vector{T}
    meass::Vector{AbstractMeasureRef}
    function MeasureExpr(cv::Vector{T}, mv::Vector{AbstractMeasureRef}) where T <: Number
        @assert length(cv) == length(mv) "Inputs need to have same length."
        @assert length(unique(mv)) == length(mv)
        @assert all(compatible.(ref, mv))
        if length(mv) > 1
            @assert all(compatible.(first(mv), mv)) "Measures do not act on same variables."
        end
        return new{T}(cv, mv)
    end
end

MP.coefficients(me::MeasureExpr) = me.coefs
measures(me::MeasureExpr) = me.meass

function Base.isequal(me1::MeasureExpr, me2::MeasureExpr)
    return all(coefficients(me1) .== coefficients(me2)) && all(measures(me1) .== measures(me2))
end

function _no_one(c::Number)
    if c == 1
        return ""
    else
        return string(c)
    end
end

function _prep_coef(c::Number)
    if c > 0
        return " + ", _no_one(abs(c))
    elseif c < 0
        return " - ", _no_one(abs(c))
    else
        return nothing
    end
end

function Base.show(io::IO, me::MeasureExpr)
    if isempty(coefficients(me))
        print(io, 0)
    else
        s, c = _prep_coef(coefficients(me)[1])
        m = measures(me)[1]
        if coefficients(me)[1] < 0
            print(io, s)
        end
        print(io, c, m)
        for (c, m) in zip(coefficients(me)[2:end], measures(me)[2:end])
            s, c = _prep_coef(c)
            print(io, s, c, m )
        end
    end
end

# All measures in a MeasureExpr act on the same variables. This is ensured by the constructor.
poly_variables(me::MeasureExpr) = poly_variables(measures(me)[1])

# compatibility

function Base.promote_rule(::Type{MeasureExpr{T}}, ::Type{<:AbstractMeasureRef}) where T
    return MeasureExpr{T}
end

function Base.convert(::Type{MeasureExpr{T}}, m::AbstractMeasureRef) where T
    return MeasureExpr([one(T)], [m])
end

function Base.promote_rule(::Type{MeasureExpr{T1}}, ::Type{MeasureExpr{T2}}) where {T1, T2}
    return MeasureExpr{promote_type(T1, T2)}
end

function Base.convert(::Type{MeasureExpr{T}}, m::MeasureExpr) where T
    return MeasureExpr(convert.(T, coefficients(m)), measures(m))
end

# Linear operations
function Base.:*(a::Number, m::AbstractMeasureRef)
    return MeasureExpr([a], [m])
end

function Base.:*(a::Number, m::MeasureExpr)
    return MeasureExpr(a.*coefficients(m), measures(m))
end

function Base.:*(m::Union{AbstractMeasureRef, MeasureExpr}, a::Number)
    return a*m
end

function Base.:/(m::Union{AbstractMeasureRef, MeasureExpr}, a::Number)
    return 1/a*m
end

function Base.:-(m::AbstractMeasureRef)
    return MeasureExpr([-1], [m])
end

function Base.:-( m::MeasureExpr)
    return MeasureExpr(-coefficients(m), measures(m))
end

function Base.sum(mev::Vector{<:MeasureExpr})
    coefs = copy(coefficients(mev[1]))
    meass = copy(measures(mev[1]))
    for me in mev[2:end]
        for (c, m) in zip(coefficients(me), measures(me))
            index = findfirst(==(m), meass)
            if isempty(index)
                push!(coefs, c)
                push!(meass, m)
            else
                coefs[index] += c
            end
        end
    end
    index = findall( x -> !(x==0.0), coefs)
    return MeasureExpr(coefs[index], meass[index])
end

function Base.sum(mev::Vector{AbstractMeasureRef})
    coefs = [1]
    meass = [mev[1]]
    for m in mev[2:end]
        index = findfirst(==(m), meass)
        if isempty(index)
            push!(coefs, 1)
            push!(meass, m)
        else
            coefs[index] += 1
        end
    end
    return MeasureExpr(coefs, meass)
end

function Base.:+(m1::Union{AbstractMeasureRef, MeasureExpr}, m2::Union{AbstractMeasureRef, MeasureExpr})
    return sum([m1, m2])
end

function Base.:-(m1::Union{AbstractMeasureRef, MeasureExpr}, m2::Union{AbstractMeasureRef, MeasureExpr})
    return sum([m1, -m2])
end

mutable struct AffineMeasureExpr{T <: Number} <: AbstractMeasureExpressionLike
    expr::MeasureExpr{T}
    cons::T
end
