function same_length(v1::AbstractVector, v2::AbstractVector)
    @assert length(v1) == length(v2) "Inputs do not have same length."
end

function same_variables(vrefs::Vector{GMPVariableRef})
    if !isempty(vrefs)
        @assert all(x -> variables(vrefs[1]) == variables(x), vrefs[2:end]) "GMPObjects do not act on same variables."
    end
end

function same_vref_type(vrefs::Vector{GMPVariableRef})
    if !isempty(vrefs)
        @assert all(x -> vref_type(vrefs[1]) == vref_type(x), vrefs[2:end]) "GMPObjects are of different type."
    end

end

"""
    ObjectExpr{C <: Number}

An ObjectExpr is a linear combination of GMPVariableRefs. As such it represends either a measure (AbstractGMPMeasure) or a continous function (AbstractGMPContinuous). 
In consequence all variables of an ObjectExpr have to refer to either AbstractGMPMeasures or AbstractGMPContinouses. In addition, these AbstractGMPObjects need to be compatible in the following sense:
  * All AbstractGMPObjects need to act on the same polynomial variables.

"""
mutable struct ObjectExpr{C <: Number} <: AbstractGMPExpr
    coefs::Vector{C}
    vars::Vector{GMPVariableRef}
    function ObjectExpr(cv::Vector{C}, mv::Vector{GMPVariableRef}) where {C <: Number}
        same_length(cv, mv)
        same_variables(mv)
        same_vref_type(mv)
        return new{C}(cv, mv)
    end
end

function ObjectExpr(c::T, v::GMPVariableRef) where {T}
    return ObjectExpr([c], [v])
end

Base.:*(a::Number, vref::GMPVariableRef) = ObjectExpr([a], [vref])
Base.:*(a::Number, e::ObjectExpr) = ObjectExpr(a.*coefficients(e), variables(e))

function Base.sum(mev::Vector{ObjectExpr{C}}) where {C}
    # TODO add early compatibiliy check: for now the sum is computed 
    # and an assertion error is generated when buildinf the GMPExpr in the return. 
    coefs = C[]
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
    index = findall( x -> x != 0.0, coefs)
    return ObjectExpr(coefs[index], vars[index])
end

function Base.sum(mev::Vector{GMPVariableRef})
    # TODO add early compatibiliy check: for now the sum is computed 
    # and an assertion error is generated when buildinf the GMPExpr in the return. 
    coefs = Int[]
    vars = GMPVariableRef[]
    for m in mev
        index = findfirst(x -> x == m, vars)
        if index isa Nothing
            push!(coefs, 1)
            push!(vars, m)
        else
            coefs[index] += c
        end
    end
    index = findall( x -> x != 0.0, coefs)
    return ObjectExpr(coefs[index], vars[index])
end

function JuMP.function_string(io, e::ObjectExpr)
    str = ""
    if length(e) == 0
        str = "0"
    else
        for (i, (c, v)) in enumerate(e)
            s, cs = _prep_coef(c)
            if (i == 1 && c > 0) 
                s = ""
            end
            str = str*s*cs*sprint(show, v)
        end 
    end
    return str
end


# compatibility

function Base.promote_rule(::Type{ObjectExpr{C}}, ::GMPVariableRef) where {C}
    return ObjectExpr{C}
end

function Base.convert(::Type{ObjectExpr{C}}, m::GMPVariableRef) where {C}
    return ObjectExpr([one(C)], [m])
end

function Base.promote_rule(::Type{ObjectExpr{C1}}, ::Type{ObjectExpr{C2}}) where {C1, C2}
    return ObjectExpr{promote_type(C1, C2)}
end

function Base.convert(::Type{ObjectExpr{C}}, m::ObjectExpr) where {C}
    return ObjectExpr(convert.(C, gmp_coefficients(m)), gmp_variables(m))
end


# The ObjectExr constructor checks that all objects in one expression have the same polynomial variables and the vref_type. 
MP.variables(e::ObjectExpr) = variables(first(gmp_variables(e)))
vref_type(e::ObjectExpr) = vref_type(first(gmp_variables(e)))

# MomentExpr

function all_vars_measures(mv::Vector{<:ObjectExpr})
    if !isempty(mv)
        @assert vref_type(first(mv)) == AbstractGMPMeasure "Polynomials can only be paired with measures"
    end
end
    
function cover_variables(c::Number, e::ObjectExpr) end
function cover_variables(c::MP.AbstractPolynomialLike, e::ObjectExpr) 
    @assert variables(c) ⊆ variables(e) "$e does not act on $c."
end

export MomentExpr
"""
    MomentExpr{T, C}

Type representing a linear combination of Moments. A Moment is the pairing of an MultivariatePolynomials.AbstractPolynomialLike and an ObjectExpr. The variables of the ObjectExpr have to be AbstractGMPMeasures which are acting on the variables of the MP.AbstractPolynomialLike.  
"""
mutable struct MomentExpr{T <: Union{Number, MP.AbstractPolynomialLike}, S <: Number} <: AbstractGMPExpr
    coefs::Vector{T}
    vars::Vector{ObjectExpr{S}}
    function MomentExpr(cv::Vector{T}, mv::Vector{ObjectExpr{S}}) where {S, T}
        same_length(cv, mv)
        all_vars_measures(mv)
        for (c, m) in zip(cv, mv)
            cover_variables(c, m)
        end
        return new{T, S}(cv, mv)
    end
end

function MomentExpr(cv::Vector{<:Union{Number, MP.AbstractPolynomialLike}}, mv::Vector{GMPVariableRef})
    return MomentExpr(cv, ObjectExpr.(ones(Int, length(mv)), mv))
end

function MomentExpr(cv::Union{Number, MP.AbstractPolynomialLike}, mv::Union{GMPVariableRef, ObjectExpr}) 
    return MomentExpr([cv], [mv])
end

degree(e::MomentExpr) = maximum(MP.maxdegree.(gmp_coefficients(e)))

Base.:*(a::Number, e::MomentExpr) = MomentExpr(a.*gmp_coefficients(e), gmp_variables(e))

function Base.sum(mev::Vector{MomentExpr{T, S}}) where {T, S}
    # TODO add early compatibiliy check: for now the sum is computed 
    # and an assertion error is generated when buildinf the GMPExpr in the return. 
    coefs = T[]
    vars = ObjectExpr{S}[]
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
    return MomentExpr(coefs[index], vars[index])
end

function JuMP.function_string(io, e::MomentExpr)
    str = ""
    if length(e) == 0
        str = "⟨0, 0⟩"
    else
        for (i, (c, v)) in enumerate(e)
            if i == 1 
                s = ""
            else 
                s = " + "
            end
            str = str*s*"⟨"*sprint(show, c)*", "*sprint(show, v)*"⟩"
        end 
    end
    return str
end

function momexp_by_measure(e::MomentExpr{T, S}) where {T, S}
    dict = Dict{AbstractGMPMeasure, promote_type(T, S)}()
    for (t, v) in e
        for (s, m) in v
            if haskey(dict, m)
                dict[m] += t*s
            else
                dict[m] = t*s
            end
        end
    end
    return dict
end

# compatibility

function Base.promote_rule(::Type{MomentExpr{T, S}}, ::Type{GMPVariableRef}) where {T, S}
    return MomentExpr{T, S}
end

function Base.convert(::Type{MomentExpr{T, S}}, m::GMPVariableRef) where {T, S, V}
    return MomentExpr([one(T)], [convert(ObjectExpr{S}, m)])
end

function Base.promote_rule(::Type{MomentExpr{T, S1}}, ::Type{ObjectExpr{S2}}) where {T, S1, S2}
    return MomentExpr{T, promote_type(S1, S2)}
end

function Base.convert(::Type{MomentExpr{T, S}}, m::ObjectExpr) where {T, S}
    return MomentExpr([one(T)], [convert(ObjectExpr{S}, m)])
end

function Base.promote_rule(::Type{MomentExpr{C1, V1}}, ::Type{MomentExpr{C2, V2}}) where {C1, C2, V1, V2}
    return MomentExpr{promote_type(T1, T2), promote_type(V1, V2)}
end

function Base.convert(::Type{MomentExpr{C, V}}, m::MomentExpr) where {C, V}
    return MomentExpr(convert.(C, gmp_coefficients(m)), convert.(ObjectExpr{V}, gmp_variables(m)))
end


export Mom
const Mom = MomentExpr
"""
    AffMomentExpr{S, T}

Type representing an affine linear combination of Moments.
"""
#=
mutable struct AffMomentExpr{T <: Union{Number, MP.AbstractPolynomialLike}, S <: Number, U <: Number} <: AbstractGMPExpressionLike
    expr::MomentExpr{T, S}
    cons::U
end
=#

const AffMomentExpr{T <: Union{Number, MP.AbstractPolynomialLike}, S <: Number, U <: Number} = 
GMPAffExpr{U, MomentExpr{T, S}}

