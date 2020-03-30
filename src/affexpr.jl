# ObjectExpr
"""
    compatible(vrefs::Vector{GMPVariableRef})

A vector of GMPVariableRefs is compatible if all objects refered to, act on the same polynomial variables and have the same var_type. 
"""
function compatible(vrefs::Vector{GMPVariableRef})
    isempty(vrefs) ? true : all(x -> variables(vrefs[1]) == variables(x), vrefs) && all( x -> x.var_type == vrefs[1].var_type, vrefs)
end

""" 
    compatible(cv::Vector{<:Number}, mv::Vector{GMPVariableRef}) 

"""
compatible(cv::Vector{<:Number}, mv::Vector{GMPVariableRef}) = length(cv) == length(mv)

"""
    ObjectExpr{C}

Type representing a linear combination of GMPVariables, i.e. mesures and continuous functions. 
Note that it is impossible to generate one linear combination with measures and continuous functions at the same time. 
"""
mutable struct ObjectExpr{C <: Number}
    coefs::Vector{C}
    vars::Vector{GMPVariableRef}
    function ObjectExpr(cv::Vector{C}, mv::Vector{GMPVariableRef})
        @assert compatible(cv, mv) "Coefficients and variables are not compatible."
        @assert compatible(mv) "Expression variables are not compatible."
        return new{C, V}(cv, mv)
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
    coefs = IntC[]
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
    return GMPExpr{C}
end

function Base.convert(::Type{ObjectExpr{C}}, m::GMPVariableRef) where {C}
    return GMPExpr([one(C)], [m])
end

function Base.promote_rule(::Type{ObjectExpr{C1}}, ::Type{ObjectExpr{C2}}) where {C1, C2}
    return GMPExpr{promote_type(C1, C2)}
end

function Base.convert(::Type{ObjectExpr{C}}, m::ObjectExpr) where {C}
    return GMPExpr(convert.(C, coefficients(m)), convert.(V, variables(m)))
end


# MomentExpr

compatible(mv::Vector{ObjectGMPExpr}) = isempty(mv) ? true : all(x -> variables(variables(first(mv))) == variables(variables(x)), vrefs) && compatible(variables.(mv))

function compatible(cv::Vector{<:MP.AbstractPolynomialLike}, mv::Vector{GMPExpr})
    if isempty(cv)
        return length(cv) == length(mv)
    else
        return length(cv) == length(mv) && variables(cv) == variables(variables(frist(mv)))
    end
end

"""
    MomentExpr{T, C}

Type representing a linear combination of Moments. A Moment is the pairing of an MultivariatePolynomials.AbstractPolynomialLike and an ObjectExpr. The variables of the ObjectExpr have to be AbstractGMPMeasures. 
"""
mutable struct MomentExpr{T <: Union{Number, MP.AbstractPolynomialLike}, S <: Number} <: AbstractGMPExpr
    coefs::Vector{T}
    vars::Vector{ObjectExpr{S}}
    function MomentExpr(cv, mv)
        @assert compatible(cv, mv) "Coefficients and variables are not compatible."
        @assert compatible(mv) "Expression variables are not compatible."
        @assert first(variables(variables(mv))).var_type isa AbstractGMPMeasure "MomentExpression need to be defined on measures."
        return new(cv, mv)
    end
end

function MomentExpr(cv::Vector{<:Union{Number, MP.AbstractPolynomialLike}}, mv::Vector{GMPVariableRef})
    return MomentExpr(cv, ObjectExpr.(ones(Int, length(mv)), mv))
end

function MomentExpr(cv::Union{Number, MP.AbstractPolynomialLike}, mv::Union{GMPVariableRef, ObejctExpr}) 
    return MomentExpr([cv], [mv])
end

Base.:*(a::Number, e::MomentExpr) = MomentExpr(a.*coefficients(e), variables(e))

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

function JuMP.function_string(io, e::ObjectExpr)
    str = ""
    if length(e) == 0
        str = "⟨0, 0⟩"
    else

        #TODO
        #
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

function Base.promote_rule(::Type{GMPExpr{C, V1}}, ::Type{V2}) where {C, V1, V2}
    return GMPExpr{C, promote_type(V1, V2)}
end

function Base.convert(::Type{GMPExpr{C, V1}}, m::V2) where {C, V1, V2}
    return GMPExpr([one(T)], [convert(V1, m)])
end

function Base.promote_rule(::Type{GMPExpr{C1, V1}}, ::Type{GMPExpr{C2, V2}}) where {C1, C2, V1, V2}
    return GMPExpr{promote_type(T1, T2), promote_type(V1, V2)}
end

function Base.convert(::Type{GMPExpr{C, V}}, m::GMPExpr) where {C, V}
    return GMPExpr(convert.(C, coefficients(m)), convert.(V, variables(m)))
end





"""
    AffMomentExpr{C}

Type representing an affine linear combination of Moments.
"""




# compatibility

function Base.promote_rule(::Type{GMPExpr{C, V1}}, ::Type{V2}) where {C, V1, V2}
    return GMPExpr{C, promote_type(V1, V2)}
end

function Base.convert(::Type{GMPExpr{C, V1}}, m::V2) where {C, V1, V2}
    return GMPExpr([one(T)], [convert(V1, m)])
end

function Base.promote_rule(::Type{GMPExpr{C1, V1}}, ::Type{GMPExpr{C2, V2}}) where {C1, C2, V1, V2}
    return GMPExpr{promote_type(T1, T2), promote_type(V1, V2)}
end

function Base.convert(::Type{GMPExpr{C, V}}, m::GMPExpr) where {C, V}
    return GMPExpr(convert.(C, coefficients(m)), convert.(V, variables(m)))
end


