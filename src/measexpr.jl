abstract type AbstractGMPScalar <: JuMP.AbstractJuMPScalar end
Base.broadcastable(t::AbstractGMPScalar) = Ref(t)


# After this only need to implement Base.:*(::Number, ::T) and Base.sum(::Vector{T}) where T <: ExprVars
const ExprVars = Union{GMPVariableRef, AbstractGMPScalar}

function Base.:*(m::ExprVars, a::Number)
    return a*m
end

function Base.:/(m::ExprVars, a::Number)
    return inv(a)*m
end

function Base.:-(m::ExprVars)
    return (-1)*m
end

function Base.:+(m1::ExprVars, m2::ExprVars)
    return sum([m1, m2])
end

function Base.:-(m1::ExprVars, m2::ExprVars)
    return sum([m1, -m2])
end

abstract type AbstractMeasureExpressionLike <: AbstractGMPScalar end

function same_length(v1::AbstractVector, v2::AbstractVector)
    @assert length(v1) == length(v2) "Inputs do not have same length."
end

function same_variables(vrefs::Vector{GMPVariableRef})
    if !isempty(vrefs)
        @assert all(x -> variables(vrefs[1]) == variables(x), vrefs[2:end]) "Measures do not act on same variables."
    end
end

export MeasExpr
"""
    MeasExpr{S}

Type representing an affine linear combination of Measures.
"""
mutable struct MeasExpr{T<:Number} <: AbstractMeasureExpressionLike
    coefs::Vector{T}
    meass::Vector{GMPVariableRef}
    function MeasExpr(coefs::Vector{T}, meass::Vector{GMPVariableRef}) where T
        same_length(coefs, meass)
        same_variables(meass)
        return new{T}(coefs, meass)
    end
end
function MeasExpr(c::T, meas::GMPVariableRef) where T
    return MeasExpr([c], [meas])
end

degree(::MeasExpr) = 0

MP.coefficients(me::MeasExpr) = me.coefs
measures(me::MeasExpr) = me.meass

MP.coefficients(me::GMPVariableRef) = [1]
measures(me::GMPVariableRef) = [me]

Base.zero(::MeasExpr{T}) where T = MeasExpr{T}(T[], GMPVariableRef[])
Base.iszero(me::MeasExpr) = length(me) == 0
Base.:(==)(me1::MeasExpr, me2::MeasExpr) = all(coefficients(me1) .== coefficients(me2)) && all(measures(me1) == measures(me2))

Base.eltype(::MeasExpr{T}) where {T} = Tuple{T,GMPVariableRef}
Base.length(me::MeasExpr) = length(coefficients(me))
Base.iterate(me::MeasExpr) = iszero(me) ? nothing : ((first(coefficients(me)), first(measures(me))), 2)
Base.iterate(me::MeasExpr, i) = length(me) < i ? nothing : ((coefficients(me)[i], measures(me)[i]), i+1)

function JuMP.function_string(mode, me::MeasExpr)
    str = ""
    if length(me) == 0
        str = "0"
    else
        for (i, (c, v)) in enumerate(me)
            if i == 1
                s = ""
                cs = string(c)
            elseif c > 0 
                s = " + "
                cs = string(c)
            else
                s = " - "
                cs = string(abs(c))
            end
            if abs(c) == 1
                cs = ""
            end
            str = str*s*cs*sprint(show, v)
        end 
    end
    return str
end

function Base.getindex(me::MeasExpr{T}, mu::GMPVariableRef) where T
    idx = findall(x -> x == mu, measures(me))
    if isempty(idx)
        return zero(T)
    else
        return sum(coefficients(me)[idx])
    end
end

Base.promote_rule(::Type{MeasExpr{T}}, ::Type{GMPVariableRef}) where T = MeasExpr{T}
Base.convert(::Type{MeasExpr{T}}, vref::GMPVariableRef) where T = MeasExpr([one(T)], [vref])
Base.promote_rule(::Type{MeasExpr{T}}, ::Type{MeasExpr{S}}) where {T, S} = MeasExpr{promote_type(T, S)}
Base.convert(::Type{MeasExpr{T}}, me::MeasExpr) where T = MeasExpr(convert.(T, coefficients(me)), measures(me))

MP.variables(me::MeasExpr) = iszero(me) ? nothing : variables(first(measures(me)))

# linear operations
Base.:*(a::Number, vref::GMPVariableRef) = MeasExpr([a], [vref])
Base.:*(a::Number, e::MeasExpr) = MeasExpr(a.*coefficients(e), measures(e))

function Base.sum(mev::Vector{GMPVariableRef})
    coefs = Int[]
    vars = GMPVariableRef[]
    for m in mev
        index = findfirst(x -> x == m, vars)
        if index isa Nothing
            push!(coefs, 1)
            push!(vars, m)
        else
            coefs[index] += 1
        end
    end
    index = findall( x -> x != 0, coefs)
    return MeasExpr(coefs[index], vars[index])
end

function Base.sum(mev::Vector{MeasExpr{T}}) where {T}
    coefs = T[]
    vars = GMPVariableRef[]
    for me in mev
        for (c, m) in me
            index = findfirst(x -> x == m, vars)
            if index isa Nothing
                push!(coefs, c)
                push!(vars, m)
            else
                coefs[index] += c
            end
        end
    end
    index = findall( x -> !iszero(x), coefs)
    return MeasExpr(coefs[index], vars[index])
end

mutable struct AffMeasExpr{S} <: AbstractMeasureExpressionLike
    expr::MeasExpr{S}
    cons::AnalyticMeasure
end

expr(ame::AffMeasExpr) = ame.expr
JuMP.constant(ame::AffMeasExpr) = ame.cons
MP.coefficients(ame::AffMeasExpr) = coefficients(expr(ame))
measures(ame::AffMeasExpr) = measures(expr(ame))

Base.:(==)(me1::AffMeasExpr, me2::AffMeasExpr) = coefficients(me1) == coefficients(me2) && measures(me1) == measures(me2) 

function JuMP.function_string(mode, ame::AffMeasExpr)
    return function_string(mode, expr(ame))*" + "*string(JuMP.constant(ame))
end

Base.:+(me::AffMeasExpr, meas::AnalyticMeasure) = AffMeasExpr(expr(me), constant(me) + meas)
Base.:+(me::MeasExpr, meas::AnalyticMeasure) = AffMeasExpr(me, meas)
Base.:+(me::GMPVariableRef, meas::AnalyticMeasure) = AffMeasExpr(convert(MeasExpr{Int}, me), meas)
Base.:-(me::Union{GMPVariableRef, MeasExpr, AffMeasExpr}, meas::AnalyticMeasure) = me + (-meas)
Base.:*(a::Number, me::AffMeasExpr) =  AffMeasExpr(a*expr(me), a*constant(me))
