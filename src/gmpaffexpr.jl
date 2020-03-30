broadcastable(t::AbstractGMPScalar) = Ref(t)

abstract type AbstractGMPExpressionLike <: AbstractGMPScalar end


"""
    abstract type AbstractGMPExpr <: AbstractGMPExpressionLike
    
Every instance of this type needs to have fields coefs and vars, or need to define the functions MP.coefficients and MP.variables. In addition the function Base.:* must be implemented for multiplication where the first element is a Number. Also Base.sum must be implemented. 
"""
abstract type AbstractGMPExpr <: AbstractGMPExpressionLike end

MP.coefficients(e::AbstractGMPExpr) = e.coefs
MP.variables(e::AbstractGMPExpr) = e.vars
Base.iszero(e::AbstractGMPExpr) = isempty(coefficients(e)) ? true : all(es -> iszero(es), coefficients(e))

function Base.getindex(e::AbstractGMPExpr, s) 
    idx = findfirst(x -> x == s, variables(e))
    idx isa Nothing ? nothing : coefficients(e)[idx]
end

Base.length(e::AbstractGMPExpr) = length(coefficients(e))
Base.iterate(e::AbstractGMPExpr) = ((first(coefficients(e)), first(variables(e))), 1)
Base.iterate(e::AbstractGMPExpr, s) = (s >= length(e)) ? nothing : ((coefficients(e)[s+1], variables(e)[s+1]), s + 1)

function Base.isequal(me1::AbstractGMPExpr, me2::AbstractGMPExpr)
    return all(coefficients(me1) .== coefficients(me2)) && all(variables(me1) .== variables(me2))
end

# Linear operations
const ExprVariables = Union{AbstractGMPVariableRef, AbstractGMPExpr}

function Base.:*(a::Number, e::ExprVariables) 
    error("Multiplication needs implementation for ($(typeof(a)), $(typeof(e))")
end

function Base.:*(m::ExprVariables, a::Number)
    return a*m
end

function Base.:/(m::ExprVariables, a::Number)
    return inv(a)*m
end

function Base.:-(m::ExprVariables)
    return (-1)*m
end

function Base.:*(::ExprVariables, ::ExpVariables)
    @error("Cannot multiply GMPExpressions. Only linear expressions are permitted.")
end
function Base.sum(mev::Vector{<:ExprVariables})
    error("Base.sum must be implemented for $(typeof(mev))")
#=
function Base.sum(mev::Vector{GMPExpr{C, V}}) where {C, V}
    # TODO add early compatibiliy check: for now the sum is computed 
    # and an assertion error is generated when buildinf the GMPExpr in the return. 
    T = eltype(mev)
    coefs = C[]
    vars = V[]
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
    index = findall( x -> x != 0.0, coefs)
    return T(coefs[index], vars[index])
end

function Base.sum(mev::Vector{S}) where S <: ExprVariables
    # TODO add early compatibiliy check: for now the sum is computed 
    # and an assertion error is generated when buildinf the GMPExpr in the return. 

    coefs = Int[]
    vars = S[]
    for m in mev
        index = findfirst(==(m), vars)
        if index isa Nothing
            push!(coefs, 1)
            push!(vars, m)
        else
            coefs[index] += 1
        end
    end
    return GMPExpr(coefs, vars)
end
=#
function Base.:+(m1::ExprVariables, m2::ExprVariables)
    return sum([m1, m2])
end

function Base.:-(m1::ExprVariables, m2::ExprVariables)
    return sum([m1, -m2])
end

i#=
mutable struct GMPAffExpr{S, T, V} <: AbstractGMPExpressionLike
    constant::S
    terms:: GMPExpr{T, V}
end

function Base.iszero(expr::GMPAffExpr)
    return iszero(expr.constant) && all(iszero, expr.terms)
end
=#
struct GMPEmptyExpr <: AbstractGMPExpressionLike end

# functions used for printing of GMPAffExpr.
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
        return "", ""
    end
end
