export GMPModel

abstract type AbstractRelaxationStatus end
struct NotRelaxed <: AbstractRelaxationStatus end
struct PrimalRelaxed <: AbstractRelaxationStatus end
struct DualRelaxed <: AbstractRelaxationStatus end
mutable struct ModelInfo
    measure_refs::Vector{JuMP.AbstractVariableRef}
    constraint_refs::Vector{AbstractGMPConstraintRef}
    variable_names::Vector{String}
    constraint_names::Vector{String}
    relax_status::AbstractRelaxationStatus
    max_degree::Int
end

ModelInfo() = ModelInfo(
                        JuMP.AbstractVariableRef[],
                        AbstractGMPConstraintRef[],
                        String[],
                        String[],
                        NotRelaxed(),
                        0)

function MOI.add_variable(mi::ModelInfo)
    push!(mi.variable_names, "")
    return length(mi.variable_names)
end

function MOI.add_constraint(mi::ModelInfo)
    push!(mi.constraint_names, "")
    return length(mi.constraint_names)
end

mutable struct ModelData
    variables::Dict{JuMP.AbstractVariableRef, MeasureVariable}
    objective::Union{Nothing, MomObj}
    constraints::Dict{AbstractGMPConstraintRef, Any}
end

ModelData() = ModelData(Dict{JuMP.AbstractVariableRef, MeasureVariable}(),
                        nothing,
                        Dict{AbstractGMPConstraintRef, Any}())


abstract type AbstractRelaxationMode end
struct PrimalFormulation <: AbstractRelaxationMode end
struct DualFormulation <: AbstractRelaxationMode end

mutable struct RelaxInfo
    mode::AbstractRelaxationMode
    relax_degree::Int
    solver::Union{Nothing, Function}
end

RelaxInfo() = RelaxInfo(DualFormulation(), 0, nothing) 

mutable struct GMPModel <: JuMP.AbstractModel
    model_info::ModelInfo
    model_data::ModelData
    relax_info::RelaxInfo
    relax_model::Union{Nothing, JuMP.Model}
end


GMPModel() = GMPModel(ModelInfo(), ModelData(), RelaxInfo(), nothing)

model_info(gmp::GMPModel) = gmp.model_info
model_data(gmp::GMPModel) = gmp.model_data
relax_info(gmp::GMPModel) = gmp.relax_info
relax_model(gmp::GMPModel) = gmp.relax_model

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
