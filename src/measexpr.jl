export MeasureExpr

const AbstractMeasureLike = Union{AbstractKnownMeasure, AbstractMeasureRef}

function compatible(meas1::AbstractMeasureLike, meas2::AbstractMeasureLike)
    if sort!(poly_variables(meas1)) == sort!(poly_variables(meas2))
        return true
    else
        return false
    end
end

mutable struct MeasureExpr{T <: Number}
    coefs::Vector{T}
    meass
    function MeasureExpr(cv::Vector{T}, mv::Vector{S}) where {T <: Number, S <: AbstractMeasureLike} 
        @assert length(cv) == length(mv) "Inputs need to have same length."
        @assert length(unique(mv)) == length(mv)
        if length(mv) > 1
            ref = first(mv)
            if !(all(compatible.(ref, mv)))
                throw(ArgumentError("Measures do not act on same variables."))
                return nothing
            end
        end
        return new{T}(cv, mv)
    end
end

Base.broadcastable(me::MeasureExpr) = Ref(me)
MP.coefficients(me::MeasureExpr) = me.coefs
measures(me::MeasureExpr) = me.meass

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
        if length(coefficients(me)) > 1
            for (c, m) in zip(coefficients(me)[2:end], measures(me)[2:end])
                s,c = _prep_coef(c)
                print(io, s, c)
                print(io, m)
            end
        end
    end
end

poly_variables(me::MeasureExpr) = poly_variables(measures(me)[1])

function Base.promote_rule(::Type{MeasureExpr{T}}, ::Type{<:AbstractMeasureLike}) where T
    return MeasureExpr{T}
end

function Base.convert(::Type{MeasureExpr{T}}, m::AbstractMeasureLike) where T
    return MeasureExpr([one(T)], [m])
end

function Base.promote_rule(::Type{MeasureExpr{T1}}, ::Type{MeasureExpr{T2}}) where {T1, T2}
    return MeasureExpr{promote_type(T1, T2)}
end

function Base.convert(::Type{MeasureExpr{T}}, m::MeasureExpr) where T
    return MeasureExpr(convert.(T, coefficients(m)), measures(m))
end

"""
Linear operations
"""

function Base.:*(a::Number, m::MeasureExpr)
    return MeasureExpr(a.*coefficients(m), measures(m))
end

function Base.:*(m::MeasureExpr, a::Number)
    return a*m
end

function Base.:/(m::MeasureExpr, a::Number)
    return 1/a*m
end

function Base.:-( m::MeasureExpr)
    return MeasureExpr(-coefficients(m), measures(m))
end

function Base.:-(m::AbstractMeasureLike)
    return MeasureExpr([-1], [m])
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

function Base.:+(m1::MeasureExpr, m2::MeasureExpr)
    return sum([m1, m2])
end

function Base.:-(m1::MeasureExpr, m2::MeasureExpr)
    return sum([m1, -m2])
end

