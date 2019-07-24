export GMPModel, add_measure!, measures, add_constraints!, add_constraint!

"""
A Generalized Moment Problem
"""
mutable struct GMPModel <:JuMP.AbstractModel
    objective::Union{Nothing, MomObj}
	constraints::Vector{MomCon}
    constraint_names::Vector{String}
	measures::Vector{Measure}

	dual::JuMP.Model
    cref::Dict{Measure, JuMP.ConstraintRef{JuMP.Model}}
    dref::Vector{JuMP.VariableRef}
    dstatus::MathOptInterface.TerminationStatusCode
end

Base.broadcastable(gmp::GMPModel) = Ref(gmp)

function GMPModel()
    model = Model()
    return GMPModel(nothing, MomCon{Int}[], String[], Measure[],
                    model, Dict{Measure,JuMP.ConstraintRef{JuMP.Model}}(), JuMP.VariableRef[], termination_status(model))
end

JuMP.object_dictionary(gmp::GMPModel) = JuMP.object_dictionary(gmp.dual)
JuMP.constraint_type(::GMPModel) = JuMP.ConstraintRef{GMPModel,Int,MomConShape}

# printing
function Base.show(io::IO, gmp::GMPModel)
    println(io, "GMPModel:")
    if gmp.objective === nothing
        println(io, "Feasibility problem:")
    else
        println(io, gmp.objective)
    end
    if !isempty(constraints(gmp))
        println(io, "s.t.")
    end
    for i in eachindex(constraints(gmp))
        con = gmp.constraints[i]
        name = gmp.constraint_names[i]
        println(io, JuMP.constraint_string(JuMP.REPLMode, name, con))
    end
    print(io, "Unknowns:")
    for μ in gmp.measures
        print(io, " $μ")
    end
    println(io)
    println(io, gmp.dstatus)

end

# add measures to GMPModel
function add_measure!(m::GMPModel, mu::Measure)
    push!(m.measures,mu)
    unique!(m.measures)
end

function add_measure!(m::GMPModel, name::String, vars::Vector{V} ;kwargs...) where V<:MP.AbstractVariable
    push!(m.measures,Measure(name,vars;kwargs...))
end

function measures(gmp::GMPModel)
    return gmp.measures
end

# add constraints to GMPModel
function constraints(gmp::GMPModel)
    return gmp.constraints
end
JuMP.name(cref::JuMP.ConstraintRef{GMPModel}) = cref.model.constraint_names[cref.index]
JuMP.constraint_object(cref::JuMP.ConstraintRef{GMPModel}) = cref.model.constraints[cref.index]

function JuMP.add_constraint(gmp::GMPModel, momcon::MomCon, name::String="noname")
    push!(constraints(gmp), momcon)
    push!(gmp.constraint_names, name)
    return ConstraintRef(gmp, length(constraints(gmp)), MomConShape())
end

function JuMP.build_constraint(_error::Function,ae::AffMomExpr,set::MOI.AbstractScalarSet)
    return MomCon(ae,set)
end

# set objective of GMPModel
function JuMP.set_objective(gmp::GMPModel, sense::MOI.OptimizationSense, mom::AbstractMomentExpressionLike)
    obj = MomObj(sense, mom)
    if isempty(setdiff(measures(obj),measures(gmp)))
        gmp.objective = obj
    else
        @error "The model does not involve $(setdiff(measures(gmp),measures(momcons)))"
    end
end
