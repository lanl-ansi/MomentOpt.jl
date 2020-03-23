# Define GMPmodel <: JuMP.AbstractModel

export GMPModel

mutable struct GMPModel <: JuMP.AbstractModel
    model_info::ModelInfo
    model_data::Dict{Symbol, Any}    
    relax_info::RelaxInfo
    relax_model::Union{Nothing, JuMP.Model}
end

GMPModel() = GMPModel(ModelInfo(), Dict{Symbol, Any}(), RelaxInfo(), nothing)

model_info(gmp::GMPModel) = gmp.model_info
JuMP.num_variables(gmp::GMPModel) = length(model_info(gmp).variable_names)
JuMP.object_dictionary(gmp::GMPModel) = gmp.model_data
relax_info(gmp::GMPModel) = gmp.relax_info
relax_model(gmp::GMPModel) = gmp.relax_model

#=
function Base.show(io::IO, gmp::GMPModel)
    println(io, "GMPModel:")
    if model_data(gmp).objective === nothing
        println(io, "Feasibility problem:")
    else
        println(io, gmp.objective)
    end
    if !isempty(model_info(gmp).constraint_refs)
        println(io, "s.t.")
    end
    for cref in model_info(gmp).constraint_refs
        con = model_data(gmp).constraints[cref]
        name = model_info(gmp).constraint_names[index(cref)]
        println(io, JuMP.constraint_string(JuMP.REPLMode, name, con))
    end
end
=#
