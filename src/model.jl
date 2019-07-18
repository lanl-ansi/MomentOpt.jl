export GMPModel, add_measure!, measures, add_constraints!, add_constraint!

"""
A Generalized Moment Problem
"""
mutable struct GMPModel <:JuMP.AbstractModel
    objective::Union{Nothing,AbstractMomentObjective}
	constraints::Vector{MomCon}
    constraint_names::Vector{String}
	measures::Vector{Measure}

	dual::JuMP.Model
    cref::Dict{Measure, JuMP.ConstraintRef{JuMP.Model}}
	dref::Dict{Any,JuMP.VariableRef}
	dstatus::MathOptInterface.TerminationStatusCode
end


function GMPModel() 
    return GMPModel(EmptyObjective(),MomCon{Int}[],String[],Measure[],Model(),Dict{Measure,JuMP.ConstraintRef{JuMP.Model}}(),Dict{Any,JuMP.VariableRef}(), termination_status(Model()))
end

JuMP.object_dictionary(gmp::GMPModel) = JuMP.object_dictionary(gmp.dual)
JuMP.constraint_type(::GMPModel) = JuMP.ConstraintRef{GMPModel,Int,MomConShape}

# printing
function Base.show(io::IO, gmp::GMPModel)
	println(io, "GMPModel:")
	println(io,gmp.objective)
	println(io, "s.t.")
	println(io, gmp.constraints)
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

function measures(gmp::GMPModel, id::Int)
	return gmp.measures[id]
end

# add constraints to GMPModel
function constraints(gmp::GMPModel)
	return gmp.constraints
end
JuMP.name(cref::JuMP.ConstraintRef{GMPModel}) = cref.model.constraint_names[cref.index]
JuMP.constraint_object(cref::JuMP.ConstraintRef{GMPModel}) = cref.model.constraints[cref.index]

function JuMP.add_constraint(gmp::GMPModel, momcon::MomCon, name::String)
    push!(constraints(gmp), momcon)
    push!(gmp.constraint_names, name)
    return ConstraintRef(gmp,length(constraints(gmp)), MomConShape())
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
