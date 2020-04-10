Base.broadcastable(t::AbstractGMPScalar) = Ref(t)

abstract type AbstractGMPExpressionLike <: AbstractGMPScalar end
struct GMPEmptyExpr <: AbstractGMPExpressionLike end

"""
    abstract type AbstractGMPExpr <: AbstractGMPExpressionLike
    
Every instance of this type needs to have fields coefs and vars, or need to define the functions gmp_coefficients and gmp_variables. In addition the function Base.:* must be implemented for multiplication where the first element is a Number. Also Base.sum must be implemented. 
"""
abstract type AbstractGMPExpr <: AbstractGMPExpressionLike end
gmp_coefficients(vref::AbstractGMPVariableRef) = [1]
gmp_coefficients(e::AbstractGMPExpr) = e.coefs
gmp_variables(vref::AbstractGMPVariableRef) = [vref]
gmp_variables(e::AbstractGMPExpr) = e.vars
Base.iszero(e::AbstractGMPExpr) = isempty(gmp_coefficients(e)) ? true : all(es -> iszero(es), gmp_coefficients(e))
JuMP.isequal_canonical(e1::AbstractGMPExpr, e2::AbstractGMPExpr) = e1 == e2

function Base.getindex(e::AbstractGMPExpr, s) 
    idx = findfirst(x -> x == s, gmp_variables(e))
    idx isa Nothing ? nothing : gmp_coefficients(e)[idx]
end

Base.length(e::AbstractGMPExpr) = length(gmp_coefficients(e))
Base.iterate(e::AbstractGMPExpr) = ((first(gmp_coefficients(e)), first(gmp_variables(e))), 1)
Base.iterate(e::AbstractGMPExpr, s) = (s >= length(e)) ? nothing : ((gmp_coefficients(e)[s+1], gmp_variables(e)[s+1]), s + 1)

function Base.:(==)(me1::AbstractGMPExpr, me2::AbstractGMPExpr)
    return all(gmp_coefficients(me1) .== gmp_coefficients(me2)) && all(gmp_variables(me1) .== gmp_variables(me2))
end

# Linear operations
const ExprVariables = Union{AbstractGMPVariableRef, AbstractGMPExpressionLike}

function Base.:*(m::ExprVariables, a::Number)
    return a*m
end

function Base.:/(m::ExprVariables, a::Number)
    return inv(a)*m
end

function Base.:-(m::ExprVariables)
    return (-1)*m
end

function Base.:+(m1::ExprVariables, m2::ExprVariables)
    return sum([m1, m2])
end

function Base.:-(m1::ExprVariables, m2::ExprVariables)
    return sum([m1, -m2])
end


mutable struct GMPAffExpr{S, V <: AbstractGMPExpr} <: AbstractGMPExpressionLike
    expr::V
    cons::S
end

expr(ae::GMPAffExpr) = ae.expr
constant(ae::GMPAffExpr) = ae.cons

function Base.:(==)(ae1::GMPAffExpr, ae2::GMPAffExpr)
    return constant(ae1) == constant(ae2) && expr(ae1) == expr(ae2)
end

function Base.:+(a::Number, v::AbstractGMPVariableRef)
    return GMPAffExpr(1*v, a)
end

function Base.:+(a::Number, e::AbstractGMPExpr)
    return GMPAffExpr(e, a)
end

function Base.:*(a::Number, ae::GMPAffExpr{S,V}) where {S, V}
    return GMPAffExpr(a*expr(ae), a*constant(ae))
end

function Base.sum(mev::Vector{GMPAffExpr{S, V}}) where {S, V}
    return GMPAffExpr(sum(expr.(mev)), sum(constant.(mev)))
end

function Base.:+(a::Number, ae::GMPAffExpr{S, V}) where {S, V}
    return GMPAffExpr(expr(ae), constant(ae)+a)
end

function Base.:+(ae::ExprVariables, a::Number) 
    return a + ae
end

function Base.:-(a::Number, ae::ExprVariables)
    return a + (-ae)
end

function Base.:-(ae::ExprVariables, a::Number)
    return ae + (-a)
end

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

function JuMP.function_string(mode, ae::GMPAffExpr)
    s, c = _prep_coef(constant(ae))
    if constant(ae) == 1
        c = string(constant(ae))
    end
    return function_string(mode, expr(ae))*s*c
end
