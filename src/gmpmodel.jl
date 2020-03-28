"""
    ModelInfo

This type contains all information to access the data of the model.
"""
mutable struct ModelInfo
    vct::Int
    cct::Int
    variable_names::Vector{AbstractString}
    constraint_names::Vector{AbstractString}
    obj_dict::Dict{Symbol, Any}
    function ModelInfo()
        new(0, 0, AbstractString[], AbstractString[], Dict{Symbol, Any}())
    end
end

"""
    ModelData

This type contains the data of the model.
"""

mutable struct ModelData
    variables::Dict{MOI.VariableIndex, AbstractGMPVariable}
    constraints::Dict{MOI.ConstraintIndex, AbstractGMPConstraint}
    objective_sense::MOI.OptimizationSense
    objective_function::GMPAffExpr
    max_degree::Int
    function ModelData()
        new(Dict{MOI.VariableIndex, AbstractGMPVariable}(),
            Dict{MOI.ConstraintIndex, AbstractGMPConstraint}(),
            MOI.FEASIBILITY_SENSE,
            GMPAffExpr{Int, AbstractGMPVariable}(), 0)
    end
end

# define types informing about the mode and status of the model.
"""
    AbstractAproximationMode

Defines whether the primal or the dual of a problem should be used for a conic approximation.
"""

export NO_APPROXIMATION_MODE, PRIMAL_APPROXIMATION_MODE, DUAL_APPROXIMATION_MODE
abstract type AbstractApproximationMode end
struct NO_APPROXIMATION_MODE <: AbstractApproximationMode end 
struct PRIMAL_APPROXIMATION_MODE <: AbstractApproximationMode end 
struct DUAL_APPROXIMATION_MODE <: AbstractApproximationMode end

export NOT_APPROXIMATED, PRIMAL_RELAXED, DUAL_STRENGTHENED, PRIMAL_STRENGTHENED, DUAL_RELAXED
abstract type AbstractApproximationStatus end
struct NOT_APPROXIMATED <: AbstractApproximationStatus end
struct PRIMAL_RELAXED <: AbstractApproximationStatus end
struct DUAL_STRENGTHENED <: AbstractApproximationStatus end
struct PRIMAL_STRENGTHENED <: AbstractApproximationStatus end
struct DUAL_RELAXED <: AbstractApproximationStatus end

mutable struct ApproximationInfo
    mode::AbstractApproximationMode
    status::AbstractApproximationStatus
    degree::Int
    solver::Union{Nothing, Function}
    function ApproximationInfo()
        new(NO_APPROXIMATION_MODE(), NOT_APPROXIMATED(), 0, nothing)
    end
end

# define GMPmodel <: JuMP.AbstractModel
# The implementation follows basically the example in JuMP/test/JuMPExtension.jl

export GMPModel

mutable struct GMPModel <: JuMP.AbstractModel
    model_info::ModelInfo
    model_data::ModelData    
    approximation_info::ApproximationInfo
    approximation_model::Union{Nothing, JuMP.Model}
    function GMPModel()
        new(ModelInfo(), ModelData(), ApproximationInfo(), nothing)
    end
end

# define internal functions for GMPModel
model_info(model::GMPModel) = model.model_info
model_data(model::GMPModel) = model.model_data
approximation_info(model::GMPModel) = model.approximation_info
approximation_model(model::GMPModel) = model.approximation_model
gmp_variables(model::GMPModel) = model_data(model).variables
gmp_constraints(model::GMPModel) = model_data(model).constraints
variable_names(model::GMPModel) = model_info(model).variable_names
constraint_names(model::GMPModel) = model_info(model).constraint_names
JuMP.object_dictionary(model::GMPModel) = model.model_info.obj_dict
JuMP.num_variables(model::GMPModel) = length(gmp_variables(model))

function JuMP.termination_status(model::GMPModel)
    if approximation_info(model).status isa NOT_APPROXIMATED
        return (NOT_APPROXIMATED(), OptimizeNodCalled)
    else
        return (approximation_info(model).status, termination_status(approximation_model(model)))
    end
end

function update_degree(m::GMPModel, degree::Int)
    if model_info(m).max_degree < degree
        model_info(m).max_degree = degree
        approximation_info(m).degree = maximum([approximation_info(m).degree, degree])
    end
end

# 
JuMP.variable_type(::GMPModel) = GMPVariableRef
JuMP.name(vref::GMPVariableRef) = variable_names(vref.model)[vref.idx]
function set_name(v::GMPVariableRef, s::String)
    owner_model(v).model_info.variable_names[v.index] = s
    return nothing
end
function JuMP.is_valid(m::GMPModel, vref::GMPVariableRef)
    return (m === vref.model &&
            vref.idx in keys(gmp_variables(m)))
end
 
function JuMP.add_variable(m::GMPModel, v::AbstractGMPVariable, name::String = "")
    model_info(m).vct += 1
    vr_type = gmp_variable_refererence_type(v)
    vref = vr_type(m, VariableIndex(model_info(m).vct))
    gmp_variables(m)[vref.index] = v
    push!(model_info(m).variable_names, name)
    return vref
end

function JuMP.delete(m::GMPModel, vref::GMPVariableRef)
    @assert JuMP.is_valid(m, vref)
    delete!(gmp_variables(m), vref.index)
    delete!(object_dictionary(m), vref.index)
end

function JuMP.variable_by_name(m::GMPModel, name::String)
    index = findall(x -> x == name, variable_names(m))
    if index isa Nothing
        return nothing
    elseif length(index) > 1
        @error "Multiple variables have the name $name."
    else
        v = object_dictionary(m)[name]
        return GMPVariableRef(m, VariableIndex(index))
    end
end

# utilities to modify variables

# define GMPConstraintRef
const GMPConstraintRef = JuMP.ConstraintRef{GMPModel, MOI.ConstraintIndex, JuMP.ScalarShape}

# define internal functions for GMPConstraintRef
JuMP.constraint_type(::GMPModel) = GMPConstraintRef

function JuMP.is_valid(m::GMPModel, cref::GMPConstraintRef)
    return (m === cref.model &&
            cref.index in keys(gmp_constraints(m)))
end

function JuMP.add_constraint(m::GMPModel, c::AbstractGMPConstraint,
                             name::String = "")
    model_info(m).cct += 1
    cref = ConstraintRef(m, MOI.ConstraintIndex(model_info(m).cct), shape(c))
    gmp_constraints(m)[cref.index] = c
    push!(constraint_names(m), name)
    update_degree(m, degree(constraint_function(c)))
    return cref
end

function JuMP.delete(m::GMPModel, cref::GMPConstraintRef)
    @assert is_valid(m, cref)
    delete!(gmp_constraints(m), cref.index)
    delete!(object_dictionary(m), cref.index)
end

function JuMP.constraint_object(cref::GMPConstraintRef)
    return gmp_constraints(cref.model)[cref.index]
end

function JuMP.num_constraints(m::GMPModel, F::Type{<:JuMP.AbstractJuMPScalar}, S::Type{<:MOI.AbstractSet})
    # TODO: is JuMP.ScalarConstraint the right thing to ask for?
    return count(con -> con isa JuMP.ScalarConstraint{F, S}, values(gmp_constraints(m)))
end

function JuMP.num_constraints(m::GMPModel, ::Type{<:Vector{F}}, S::Type{<:MOI.AbstractSet}) where F<:JuMP.AbstractJuMPScalar
    # TODO: is JuMP.VectorConstraint the right thing to ask for?
    return count(con -> con isa JuMP.VectorConstraint{F, S}, values(gmp_constraints(m)))
end

JuMP.name(cref::GMPConstraintRef) = constraint_names(cref.model)[cref.index]

function JuMP.set_name(cref::GMPConstraintRef, name::String)
    constraint_names(cref.model)[cref.index] = name
end

function JuMP.constraint_by_name(m::GMPModel, name::String)
    index = findall(x -> x == name, constraint_names(m))
    if index isa Nothing
        return nothing
    elseif length(index) > 1
        @error "Multiple constraints have the name $name."
    else
        return GMPConstraintRef(m, MOI.ConstraintIndex(index))
    end
end

# Objective
function JuMP.set_objective(m::GMPModel, sense::MOI.OptimizationSense,
                            f::GMPAffExpr)
    model_data(m).objective_sense = sense
    model_data(m).objective_function = f
    update_degree(m, degree(f))    
end

JuMP.objective_sense(m::GMPModel) = model_data(m).objective_sense

function JuMP.set_objective_sense(m::GMPModel, sense::MOI.OptimizationSense)
    model_data(m).objective_sense = sense
end

JuMP.objective_function(m::GMPModel) = model_data(m).objective_function
JuMP.objective_function_type(m::GMPModel) = typeof(objective_function(m))

function JuMP.objective_function(m::GMPModel, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(objective_function(m) isa FT)
        throw(InexactError(:objective_function, FT, objective_function_type(m)))
    end
    return model_data(m).objective_function::FT
end

# Show
function JuMP.show_backend_summary(io::IO, model::GMPModel)
    println(io)
    println(io, "Approxmation mode: ", approximation_info(model).mode)
    println(io, "Approximation status: ", approximation_info(model).status)
    println(io, "Maximum degree of data: ", model_data(model).max_degree)
    println(io, "Degree for approximation ", approximation_info(model).degree)
    # The last print shouldn't have a new line
    if approximation_info(model).solver isa Nothing
        print(io, "Solver for approximation: ")
    else
        print(io, "Solver for approximation: ", approximation_info(model).solver)
    end
end

function JuMP.show_objective_function_summary(io::IO, m::GMPModel)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(m))
end

function JuMP.objective_function_string(print_mode, m::GMPModel)
    return JuMP.function_string(print_mode, objective_function(m))
end

_plural(n) = (isone(n) ? "" : "s")

function JuMP.show_constraints_summary(io::IO, m::GMPModel)
    n = length(gmp_constraints(m))
    print(io, "Constraint", _plural(n), ": ", n)
end

function JuMP.constraints_string(print_mode, m::GMPModel)
    strings = String[]
    # Sort by creation order, i.e. MOI.ConstraintIndex value
    constraints = sort(collect(gmp_constraints(m)), by = c -> c.first.value)
    for (index, constraint) in constraints
        push!(strings, JuMP.constraint_string(print_mode, constraint))
    end
    return strings
end
