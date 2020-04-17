# ApproximationSequences connect GMPObjects to JuMPvariables
export ApproximationSequence
"""
    ApproximationSequence{T <: AbstractGMPObject}

An approximation sequence stores necessary data to build the approximation of the object.
"""
mutable struct ApproximationSequence{T <: AbstractGMPObject}
    model::JuMP.Model
    object::T
    dict::Dict{MP.AbstractPolynomialLike, Any}
    function ApproximationSequence(o::T; model = JuMP.Model()) where T <: AbstractGMPObject
        new{T}(model, o, Dict{MP.AbstractPolynomialLike, Any}())
    end
end

Base.broadcastable(ms::ApproximationSequence) = Ref(ms)

function gmp_object(ms::ApproximationSequence)
    return ms.object
end

function Base.show(io::IO, ms::ApproximationSequence)
    println(io, "ApproximationSequence:")
    for (k, v) in ms.dict
        println("$k => $v")
    end
end

Base.length(ms::ApproximationSequence) = length(ms. dict)

function Base.getindex(ms::ApproximationSequence, m::MP.AbstractPolynomialLike)
    if !haskey(ms.dict, m) 
        if MOI.get(gmp_object(ms), ApproximationType()) isa EXACT_APPROXIMATION
            ms.dict[m] = integrate(m, gmp_object(ms))
        else
            cc, mm = MB.change_basis(m, MOI.get(gmp_object(ms), ApproximationBasis()))
            if length(cc) == 1 && mm[1] == m
                ms.dict[m] = @variable ms.model
            else
                ms.dict[m] = 0
                for i in 1:length(cc)
                    haskey(ms.dict, mm[i]) ? nothing : ms.dict[mm[i]] = @variable ms.model 
                    ms.dict[m] += cc[i] * ms.dict[mm[i]]
                end
            end
        end
    end
    return ms.dict[m]
end

function Base.getindex(ms::ApproximationSequence, m::Number)
    return getindex(ms::ApproximationSequence, m*_mono_one(MOI.get(gmp_object(ms), Variables())))
end

Base.getindex(ms::ApproximationSequence, v::AbstractArray) = getindex.(ms, v)

function approximation_sequence(model::JuMP.Model, o::AbstractGMPObject, degree::Int)
    ms = ApproximationSequence(o; model = model)
    basis = maxdegree_basis(o, degree) #TODO change for sparsity or symmetry
    ms[MP.monomials(basis)] #ApproximationSequences are generated lazily.
    return ms
end

export GMPObjectApproximation
"""
    GMPObjectApproximation

Contains approximation and references to constraints, justifying the approximation.
"""
mutable struct GMPObjectApproximation
    approximation::ApproximationSequence
    justification::Vector{JuMP.ConstraintRef}
end

function Base.show(io::IO, ap::GMPObjectApproximation)
    print(io, "Approximation with $(length(ap.justification)) crefs")
end

function init_approximation(model::GMPModel, degree::Int)
    m = approximation_model(model)
    for (i, o) in gmp_variables(model)
        approx_vrefs(model)[i] = PrimalDualRef(GMPObjectApproximation(approximation_sequence(m, o.v, degree), GMPConstraintRef[]), nothing)
    end
    return nothing
end

function approximation_sequence(model::GMPModel, index::Int)
    return approx_vrefs(model)[index].primal.approx
end

"""
    riesz(::GMPModel, ::GMPVariableRef, p)
    riesz(::GMPModel, ::(Aff)MomentExpr)

Evaluates a moment expression with respect to the corresponding moment sequences.
"""
function riesz(m::GMPModel, index::Int, p::MP.AbstractPolynomialLike)
    cc, mm = MB.change_basis(p, MOI.get(gmp_variables(m)[index].v, ApproximationBasis()))
    ms = approximation_sequence(m, index)
    return sum( cc[i]*ms[mm[i]] for i in 1:length(cc))
end

function riesz(ms::GMPModel, index::Int, m::Number)
    return riesz(ms, index, m*_mono_one(MOI.get(gmp_variables(ms)[index].v, Variables())))
end

function riesz(model::GMPModel, me::MomentExpr)
    return sum(c*riesz(model, index(m), cv) for (cv, mv) in me for (c, m) in mv)
end

riesz(ms::GMPModel, me::AffMomentExpr) = constant(me) + riesz(expr(me))

