broadcastable(t::AbstractGMPScalar) = Ref(t)

abstract type AbstractGMPExpressionLike <: AbstractGMPScalar end
const ExprVariables = Union{AbstractGMPVariableRef, AbstractGMPExpressionLike}

mutable struct GMPExpr{C <: Number, V <: ExprVariables} <: AbstractGMPExpressionLike
    coefs::Vector{C}
    vars::Vector{V}
    function GMPExpr(cv::Vector{C}, mv::Vector{V}) where {C <: Number, V <: ExprVariables}
        @assert length(cv) == length(mv) "Inputs need to have same length."
        @assert length(unique(mv)) == length(mv)
        @assert all(compatible.(ref, mv))
        if length(mv) > 1
            @assert all(compatible.(first(mv), mv)) "Measures do not act on same variables."
        end
        return new{C, V}(cv, mv)
    end
end

MP.coefficients(e::GMPExpr) = e.coefs
MP.variables(e::GMPExpr) = e.vars
Base.iszero(e::GMPExpr) = isempty(coefficients(e)) ? true : all(es -> iszero(es), coefficients(e))

function Base.getindex(e::GMPExpr{C, V}, s::V) where {C, V} 
    idx = findfirst(x -> x == s, variables(e))
    idx isa Nothing ? nothing : coefficients(e)[idx]
end
Base.iterate(e::GMPExpr) = ((first(coefficients(e), first(variables(e))), 1)
Base.iterate(e::GMPExpr, s) = (coefficients(e)[s+1], variables(e)[s+1], s + 1)
Base.length(e::GMPExpr) = length(coefficientes(e))
JuMP.variable_ref_type(::GMPExpr{C, V}) where {C, V} = V

function Base.isequal(me1::GMPExpr, me2::GMPExpr)
    return all(coefficients(me1) .== coefficients(me2)) && all(variables(me1) .== variables(me2))
end

# compatibility

function Base.promote_rule(::Type{GMPExpr{C, V1}}, ::Type{V2}) where {C, V1, V2}
    return GMPExpr{C, promote_type(V1, V2)}
end

function Base.convert(::Type{GMPExpr{C, V}}, m::V) where {C, V}
    return GMPExpr([one(T)], [m])
end

function Base.promote_rule(::Type{GMPExpr{C1, V1}}, ::Type{GMPExpr{C2, V2}}) where {C1, C2, V1, V2}
    return GMPExpr{promote_type(T1, T2), promote_type(V1, V2)}
end

function Base.convert(::Type{GMPExpr{C, V}}, m::GMPExpr) where {C, V}
    return GMPExpr(convert.(C, coefficients(m)), convert.(V, variables(m)))
end

# Linear operations

function Base.:*(a::Number, m::ExprVariables)
    # The function will probabably infer the wrong type dependingon S <: ExprVariables.
    # Therefore multiplications should be defined for all implementations of GMPExpr.    
    return GMPExpr([a], [m])
end

function Base.:*(a::S, m::GMPExpr{T, V}) where {S <: Number, T <: Number, V <: ExprVariables}
    return GMPExpr{promote_type(S,T), V}(a.*coefficients(m), variables(m))
end

function Base.:*(m::Union{ExprVariables, GMPExpr}, a::Number)
    return a*m
end

function Base.:/(m::Union{ExprVariables, GMPExpr}, a::Number)
    return inv(a)*m
end

function Base.:-(m::GMPExpr)
    return (-1)*m
end

function Base.sum(mev::Vector{GMPExpr{C, V}}) where {C, V}
    coefs = C[]
    vars = V[]
    for me in mev
        for (c, m) in me
            index = findfirst(==(m), vars)
            if isempty(index)
                push!(coefs, c)
                push!(vars, m)
            else
                coefs[index] += c
            end
        end
    end
    index = findall( x -> !(x != 0.0), coefs)
    return GMPExpr{C, V}(coefs[index], vars[index])
end

function Base.sum(mev::Vector{S}) where S <: ExprVariables
    coefs = Int[]
    vars = S[]
    for m in mev
        index = findfirst(==(m), vars)
        if isempty(index)
            push!(coefs, 1)
            push!(vas, m)
        else
            coefs[index] += 1
        end
    end
    return GMPExpr{Int, S}(coefs, meass)
end

function Base.:+(m1::Union{ExprVariables, GMPExpr}, m2::Union{ExprVariables, GMPExpr})
    return sum([m1, m2])
end

function Base.:-(m1::Union{ExprVariables, GMPExpr}, m2::Union{ExprVariables, GMPExpr})
    return sum([m1, -m2])
end


mutable struct GMPAffExpr{S, T, V} <: AbstractGMPExpressionLIke
    constant::S
    terms:: GMPExpr{T, V}
end

function Base.iszero(expr::GMPAffExpr)
    return iszero(expr.constant) && all(iszero, expr.terms)
end
