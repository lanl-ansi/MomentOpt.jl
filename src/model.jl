export GMPModel, add_measure!, measures, add_constraints!, add_constraint!

"""
A Generalized Moment Problem
"""
mutable struct GMPModel <:JuMP.AbstractModel
    objective::Union{Nothing,AbstractMomentObjective}
	constraints::Vector{<:AbstractMomentConstraint}
	measures::Vector{Measure}

	dual::JuMP.Model
    cref::Dict{Measure, JuMP.ConstraintRef{JuMP.Model}}
	dref::Dict{Any,JuMP.VariableRef}
	dstatus::MathOptInterface.TerminationStatusCode
end


function GMPModel() 
    return GMPModel(EmptyObjective(),MomCon{Int,Int}[],Measure[],Model(),Dict{Measure,JuMP.ConstraintRef{JuMP.Model}}(),Dict{Any,JuMP.VariableRef}(), termination_status(Model()))
end


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
	append!(m.measures,[mu])
	unique!(m.measures)
end

function add_measure!(m::GMPModel, name::String, vars::Vector{V} ;kwargs...) where V<:MP.AbstractVariable
	append!(m.measures,[Measure(name,vars;kwargs...)])
end

function measures(gmp::GMPModel)
	return gmp.measures
end

function measures(gmp::GMPModel, id::Int)
	return gmp.measures[id]
end

# add constraints to GMPModel
function add_constraints!(gmp::GMPModel, momcons::Vector{M}) where M<:AbstractMomentConstraint
		if length(gmp.constraints)==0
			gmp.constraints = momcons
		else
			try convert(typeof(gmp.constraints),momcons)
			catch error 
			gmp.constraints = convert(Vector{M},gmp.constraints)
			end
			append!(gmp.constraints, momcons)			
		end
end

function add_constraint!(gmp::GMPModel,momcon::M) where M<:AbstractMomentConstraint
	add_constraints!(gmp,[momcon])
end

function constraints(gmp::GMPModel)
	return gmp.constraints
end

function JuMP.set_objective(gmp::GMPModel,sense::MOI.OptimizationSense,mom::AbstractMomentExpressionLike)
    obj = MomObj(sense, mom)
    if isempty(setdiff(measures(obj),measures(gmp)))
        gmp.objective = obj
    else
        @error "The model does not involve $(setdiff(measures(gmp),measures(momcons)))"
    end
end
