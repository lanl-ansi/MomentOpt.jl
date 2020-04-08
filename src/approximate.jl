export approximate!
"""
    approximate!(m::GMPModel)

Computes an approximation of m.
"""
function approximate!(model::GMPModel)
    if approximation_info(model).mode isa NO_APPROXIMATION_MODE
        @warn("No approximation mode has been set.") #TODO remove default to Primal
    else
        generate_approximation_model(model)
        set_optimizer(approximation_model(model), approximation_info(model).solver)
        optimize!(approximation_model(model))
        #TODO maybe do some postprocessing
    end
    return model
end


# Moment sequences connect GMPObjects to JuMPvariables
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
    approx::ApproximationSequence
    crefs::Vector{JuMP.ConstraintRef}
end

function Base.show(io::IO, ap::GMPObjectApproximation)
    print(io, "Approximation with $(length(ap.crefs)) crefs")
end

function init_approximation(model::GMPModel, degree::Int)
    m = approximation_model(model)
    for (i, o) in gmp_variables(model)
        approximation_info(model).approx_vrefs[i] = GMPObjectApproximation(approximation_sequence(m, o.v, degree), GMPConstraintRef[])
    end
    return nothing
end

function approximation_sequence(model::GMPModel, index::Int)
    return approximation_info(model).approx_vrefs[index].approx
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

# inner approximations of positive polynomials

function QuadraticModuleMonomials(o::AbstractGMPObject, degree::Int)
    eqs = equalities(MOI.get(o, BSASet()))
    ineqs = [_mono_one(MOI.get(o, Variables())), inequalities(MOI.get(o, BSASet()))...]
    QM = Dict(:eq => Dict(), :ineq => Dict())
    # TODO reflect basis of Approximation
    for eq in eqs
        QM[:eq][eq] = monomials(MOI.get(o, Variables()), 0:degree-maxdegree(eq))
    end
    for ineq in ineqs 
        QM[:ineq][ineq] = monomials(MOI.get(o, Variables()), 0:Int(floor((degree-maxdegree(ineq))/2)))
    end
    return QM
end

function MomentPutinar(model::GMPModel, degree::Int)
    for (i, v_approx) in approximation_info(model).approx_vrefs
        if MOI.get(gmp_variables(model)[i].v, ApproximationType()) isa MomentOpt.DefaultApproximation
            qm = QuadraticModuleMonomials(gmp_variables(model)[i].v, degree)
            for (eq, mults) in qm[:eq]
                for mult in mults
                    cref =  @constraint approximation_model(model) riesz(model, i, eq*mult) == 0
                    push!(v_approx.cref, cref)
                end
            end
            for (ineq, mult) in qm[:ineq]
                if length(mult) == 1
                    moiset = MOI.GreaterThan(0.0)
                    mult = first(mult)
                else
                    moiset = PSDCone()
                end
                cref = @constraint approximation_model(model) riesz.(model, i, ineq*mult*transpose(mult)) in moiset
                push!(v_approx.crefs, cref)
            end
        end
    end
    return nothing
end


function generate_approximation_model(model::GMPModel)
    return measure_relaxation_model(model::GMPModel) # TODO be general
end

function measure_relaxation_model(model::GMPModel)
    model.approximation_model = Model()
    degree = approximation_info(model).degree
    # initiate moments
    init_approximation(model, degree)
    # add substitutions
    # add measure condition to moments
    MomentPutinar(model, degree)

    # from now on no new variables shoud be added. 

    # add objective
    @objective approximation_model(model) objective_sense(model) riesz(model, objective_function(model)) 
    # add moment constraints
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MomentConstraintShape
            approximation_info(model).approx_crefs[i] = @constraint approximation_model(model) riesz(model, jump_function(con)) in moi_set(con)
        end
    end
    # add measure constraints
    return nothing
end

