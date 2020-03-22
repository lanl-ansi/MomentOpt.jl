# Define abstract types and backend for GMPModel.

abstract type AbstractMeasureRef <: JuMP.AbstractVariableRef end

abstract type AbstractRelaxationStatus end
struct NotRelaxed <: AbstractRelaxationStatus end
struct PrimalRelaxed <: AbstractRelaxationStatus end
struct DualRelaxed <: AbstractRelaxationStatus end

mutable struct ModelInfo
    measure_refs::Vector{AbstractMeasureRef}
    constraint_refs::Vector{JuMP.ConstraintRef}
    variable_names::Vector{String}
    constraint_names::Vector{String}
    relax_status::AbstractRelaxationStatus
    max_degree::Int
end

ModelInfo() = ModelInfo(
                        AbstractMeasureRef[],
                        JuMP.ConstraintRef[],
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

abstract type AbstractRelaxationMode end
struct PrimalFormulation <: AbstractRelaxationMode end
struct DualFormulation <: AbstractRelaxationMode end

mutable struct RelaxInfo
    mode::AbstractRelaxationMode
    relax_degree::Int
    solver::Union{Nothing, Function}
end

RelaxInfo() = RelaxInfo(DualFormulation(), 0, nothing) 
