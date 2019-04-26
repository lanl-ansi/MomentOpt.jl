export GMPModel, add_measure!, measures, add_constraints!, add_constraint!

"""
GMP 
"""
mutable struct GMPModel
	objective::OBJ where OBJ<:AbstractMomentObjective
	constraints::Vector{<:AbstractMomentConstraint}
	measures::Vector{Measure}

	dual::JuMP.Model		
	cref::Dict{Measure,Any}
	dref::Dict{Any,Any}
	dstatus::MathOptInterface.TerminationStatusCode
end


function GMPModel() 
	return GMPModel(EmptyObjective(),CMomCon{Int,Int}[],Measure[],Model(),Dict{Measure,Any}(),Dict{Any,Any}(), termination_status(Model()))
end


# printing
function Base.show(io::IO,gmp::GMPModel)
	println("GMPModel:")
	println(io,gmp.objective)
	println("s.t.")
	println(io,gmp.constraints)
	print("Unknowns: ")
	for μ in gmp.measures
	print("$μ ")
	end
	println()
	println(io,gmp.dstatus)
	
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


function add_objective!(gmp::GMPModel, obj::MO) where MO<:AbstractMomentObjective
	if isempty(setdiff(measures(obj),measures(gmp)))
	gmp.objective = obj
	else
		println("GMP does not involve $(setdiff(measures(gmp),measures(momcons)))")
	end
end

function add_objective!(gmp::GMPModel,sense::Symbol,mom::M) where M<:AbstractMomentExpressionLike
	obj = MomObj(sense, mom)
	add_objective!(gmp,obj)
end

function add_objective!(gmp::GMPModel,sense::Symbol,pol::PT,mu::Measure) where {M<:AbstractMomentExpressionLike, PT<:AbstractPolynomialLike}
	obj = MomObj(sense, Mom(pol,mu))
	add_objective!(gmp,obj)
end

function add_objective!(gmp::GMPModel,sense::Symbol,mu::Measure,pol::PT) where {M<:AbstractMomentExpressionLike, PT<:AbstractPolynomialLike}
	obj = MomObj(sense, Mom(pol,mu))
	add_objective!(gmp,obj)
end
	
function add_objective!(gmp::GMPModel,sense::Symbol,pol::T,mu::Measure) where {M<:AbstractMomentExpressionLike, T<:Number}
	obj = MomObj(sense, CMom(pol,mu))
	add_objective!(gmp,obj)
end

function add_objective!(gmp::GMPModel,sense::Symbol,mu::Measure,pol::T) where {M<:AbstractMomentExpressionLike, T<:Number}
	obj = MomObj(sense, CMom(pol,mu))
	add_objective!(gmp,obj)
end
	
