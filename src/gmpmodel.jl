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
    variables::Dict{Int, AbstractGMPVariable}
    constraints::Dict{Int, AbstractGMPConstraint}
    objective_sense::MOI.OptimizationSense
    objective_function::AbstractMomentExpressionLike
    max_degree::Int
    function ModelData()
        new(Dict{Int, AbstractGMPVariable}(),
            Dict{Int, AbstractGMPConstraint}(),
            MOI.FEASIBILITY_SENSE,
            zero(MomentExpr{Int}), 0)
    end
end

# define types informing about the mode and status of the model.
"""
    AbstractAproximationMode

Defines whether the primal or the dual of a problem should be used for a conic approximation.
"""
abstract type AbstractApproximationMode end

export PRIMAL_RELAXATION_MODE, DUAL_STRENGTHEN_MODE

abstract type AbstractPrimalMode <: AbstractApproximationMode end
abstract type AbstractDualMode <: AbstractApproximationMode end
struct PRIMAL_RELAXATION_MODE <: AbstractPrimalMode end 
struct DUAL_STRENGTHEN_MODE <: AbstractDualMode end

struct VarApprox{S, T}
    value::MM.Measure
    p_just::S
    d_just::T
end

struct ConApprox{S, T}
    p_res::S
    dual::T
end

mutable struct ApproximationInfo
    mode::AbstractApproximationMode
    degree::Int
    approx_vrefs::Dict{Int, VarApprox}
    approx_crefs::Dict{Int, ConApprox}
    function ApproximationInfo()
        new(DUAL_STRENGTHEN_MODE(), 0, Dict{Int, VarApprox}(), Dict{Int, ConApprox}())
    end
end

# define GMPmodel <: JuMP.AbstractModel
# The implementation follows basically the example in JuMP/test/JuMPExtension.jl
export GMPModel
"""
    GMPModel

This JuMP.AbstractModel is able to model linear obtimization problems where the variables are abstract measures or continuous functions.
"""
mutable struct GMPModel <: AbstractGMPModel
    model_info::ModelInfo
    model_data::ModelData    
    approximation_info::ApproximationInfo
    approximation_model::JuMP.Model
    function GMPModel()
        new(ModelInfo(), ModelData(), ApproximationInfo(), Model())
    end
end
Base.broadcastable(model::GMPModel) = Ref(model)
# define internal functions for GMPModel
model_info(model::GMPModel) = model.model_info
model_data(model::GMPModel) = model.model_data
approximation_info(model::GMPModel) = model.approximation_info
export approximation_mode
approximation_mode(model::GMPModel) = approximation_info(model).mode
approx_vrefs(model::GMPModel) = approximation_info(model).approx_vrefs
approx_crefs(model::GMPModel) = approximation_info(model).approx_crefs
export approximation_model
approximation_model(model::GMPModel) = model.approximation_model
gmp_variables(model::GMPModel) = model_data(model).variables
gmp_constraints(model::GMPModel) = model_data(model).constraints
variable_names(model::GMPModel) = model_info(model).variable_names
constraint_names(model::GMPModel) = model_info(model).constraint_names
model_degree(model::GMPModel) = model_data(model).max_degree
export approximation_degree
approximation_degree(model::GMPModel) = approximation_info(model).degree

function JuMP.set_optimizer(model::GMPModel, optimizer)
    set_optimizer(approximation_model(model), optimizer)
end

function GMPModel(optimizer_factory)
    model = GMPModel()
    set_optimizer(model, optimizer_factory)
    return model
end

JuMP.solver_name(model::GMPModel) = solver_name(approximation_model(model))
JuMP.num_variables(model::GMPModel) = length(gmp_variables(model))
JuMP.num_nl_constraints(::GMPModel) = 0
JuMP.object_dictionary(model::GMPModel) = model_info(model).obj_dict
JuMP.raw_status(model::GMPModel) = raw_status(approximation_model(model))
JuMP.set_optimize_hook(model::GMPModel, f) = (approximation_model(model).optimize_hook = f)
JuMP.solve_time(model::GMPModel) = solve_time(approximation_model(model))

function JuMP.set_optimizer_attribute(model::GMPModel, args...)
    set_optimizer_attribute(approximation_model(model), args...)
end
function JuMP.set_optimizer_attributes(model::GMPModel, args...)
    set_optimizer_attributes(approximation_model(model), args...)
end

function JuMP.get_optimizer_attribute(model::GMPModel, name::String)
    return get_optimizer_attribute(approximation_model(model), MOI.RawParameter(name))
end

function JuMP.set_silent(model::GMPModel)
    return MOI.set(approximation_model(model), MOI.Silent(), true)
end

function JuMP.unset_silent(model::GMPModel)
    return MOI.set(approximation_model(model), MOI.Silent(), false)
end

JuMP.time_limit_sec(model::GMPModel) = time_limit_sec(approximation_model(model))

function update_degree(m::GMPModel, degree::Int)
    if model_data(m).max_degree < degree
        model_data(m).max_degree = degree
        approximation_info(m).degree = maximum([approximation_info(m).degree, 2*ceil(degree/2)])
    end
end

export set_approximation_degree
"""
    set_approximation_degree(model::GMPModel, degree::Int)

Manually set the degree of the approximation to a custom value. The value cannot be less than the degree of the data. 
"""
function set_approximation_degree(model::GMPModel, degree::Int)
    if approximation_info(model).degree > degree && model_data(model).max_degree > degree
        approximation_info(model).degree = minimum([approximation_info(model).degree, 2*ceil(model_data(model).max_degree/2)])

        @warn "Requested approximation degree $degree is too low to cover all data. The approximation degree has been set to the minimal value possible and now is $(approximation_info(model).degree)."
    else
        approximation_info(model).degree = maximum([approximation_info(model).degree, 2*ceil(degree/2)])
    end
    return nothing
end

export set_approximation_mode
"""
    set_approximation_mode(m::GMPModel, mode::AbstractApproximationMode) 

Set the approximation mode of m to mode. Depending on the solver either the primal or the dual formulation might be more efficient. It might also be interesting to change the mode when numerical issues are encountered.
Possible values for mode are:
  * PRIMAL_RELAXATION_MODE()
  * DUAL_STRENGTHEN_MODE() (default)
  """
function set_approximation_mode(m::GMPModel, mode::AbstractApproximationMode) 
    approximation_info(m).mode = mode
    return m
end

# referencing variables
function GMPVariableRef(m::GMPModel, v::AbstractGMPVariable)
    model_info(m).vct += 1
    return GMPVariableRef(m, model_info(m).vct)
end

JuMP.name(vref::GMPVariableRef) = variable_names(owner_model(vref))[index(vref)]

function JuMP.set_name(v::GMPVariableRef, s::String)
    idx = findall(x -> x == JuMP.name(v), variable_names(owner_model(v)))
    if length(idx) == 1
        delete!(object_dictionary(owner_model(v)), Symbol(JuMP.name(v)))
    end
    variable_names(owner_model(v))[index(v)] = s
    return nothing
end

function JuMP.variable_by_name(m::GMPModel, name::String)
    idx = findall(x -> x == name, variable_names(m))
    if isempty(idx)
        return nothing
    elseif length(idx) > 1
        throw(ErrorException("Multiple variables have the name $name."))
    else
        return object_dictionary(m)[Symbol(name)]
    end
end

function JuMP.is_valid(m::GMPModel, vref::GMPVariableRef)
    return (m === owner_model(vref) &&  index(vref) in keys(gmp_variables(m)))
end

function JuMP.add_variable(m::GMPModel, v::AbstractGMPVariable, name::String = "")
    vref = GMPVariableRef(m, v)
    gmp_variables(m)[index(vref)] = v
    push!(variable_names(m), name)
    update_degree(m, maxdegree(bsa_set(v.v)))
    return vref
end

function JuMP.delete(m::GMPModel, vref::GMPVariableRef)
    @assert JuMP.is_valid(m, vref)
    delete!(gmp_variables(m), index(vref))
    idx = findall(x -> x == JuMP.name(vref), variable_names(m))
    if length(idx) == 1
        delete!(object_dictionary(m), Symbol(JuMP.name(vref)))
    end
    variable_names(m)[index(vref)] = ""
    return nothing
end


function vref_object(vref::GMPVariableRef)
    return gmp_variables(owner_model(vref))[vref.index].v
end

function MP.variables(vref::GMPVariableRef)
    return variables(vref_object(vref))
end

export support
"""
    support(vref::GMPVariableRef)

Retruns the support a measure is defined on. 
"""
function support(vref::GMPVariableRef)
    return bsa_set(vref_object(vref))
end

export domain
"""
    domain(vref::GMPVariableRef)

Retruns the domain a continuous is defined on. 
"""
function domain(vref::GMPVariableRef)
    return bsa_set(vref_object(vref))
end

export approx_basis
"""
    approx_basis(vref::GMPVariableRef)

Returns the basis used to approximate vref.
"""
function approx_basis(vref::GMPVariableRef)
    return approx_basis(vref_object(vref))
end

# referencing constraints 
function JuMP.is_valid(m::GMPModel, cref::GMPConstraintRef)
    return (m === cref.model &&
            cref.index in keys(gmp_constraints(m)))
end

function JuMP.constraint_string(print_mode, ref::GMPConstraintRef; in_math_mode = false)
    return JuMP.constraint_string(print_mode, JuMP.name(ref), constraint_object(ref), in_math_mode = in_math_mode)
end

function Base.show(io::IO, ref::GMPConstraintRef)
    print(io, constraint_string(REPLMode, ref))
end

function JuMP.constraint_object(cref::GMPConstraintRef)
    return gmp_constraints(cref.model)[cref.index]
end

function JuMP.jump_function(cref::GMPConstraintRef)
    return jump_function(constraint_object(cref))
end

function JuMP.moi_set(cref::GMPConstraintRef)
    return moi_set(constraint_object(cref))
end

function JuMP.add_constraint(m::GMPModel, c::AbstractGMPConstraint,
                             name::String = "")
    model_info(m).cct += 1
    cref = GMPConstraintRef(m, model_info(m).cct, shape(c))
    gmp_constraints(m)[cref.index] = c
    push!(constraint_names(m), name)
    update_degree(m, degree(jump_function(c)))
    return cref
end

function JuMP.add_constraint(m::GMPModel, c::ScalarConstraint{<:Union{AffMomentExpr, MomentExpr}, <:MOI.AbstractScalarSet}, name::String = "") where {S}
    add_constraint(m, MomentConstraint(jump_function(c), moi_set(c)), name)
end

function JuMP.add_constraint(m::GMPModel, c::ScalarConstraint{<:Union{AffMeasExpr, MeasExpr}, <:MOI.AbstractScalarSet}, name::String = "")
    @assert moi_set(c) isa MOI.EqualTo "Only support equality measure constraints."
    @assert MOI.constant(moi_set(c)) == 0 "A measure cannot be equal to a number. "
    return add_constraint(m, MeasureConstraint(jump_function(c)), name)
end

function JuMP.delete(m::GMPModel, cref::GMPConstraintRef)
    @assert is_valid(m, cref)
    delete!(gmp_constraints(m), index(cref))
    idx = findall(x -> x == JuMP.name(cref), constraint_names(m))
    if length(idx) == 1
        delete!(object_dictionary(m), Symbol(JuMP.name(cref)))
    end
    constraint_names(m)[index(cref)] = ""
    return nothing
end

function JuMP.num_constraints(m::GMPModel, T::Type{<:AbstractGMPConstraint})
    return count(con -> con isa T, values(gmp_constraints(m)))
end

JuMP.name(cref::GMPConstraintRef) = constraint_names(cref.model)[cref.index]

function JuMP.set_name(cref::GMPConstraintRef, name::String)
    idx = findall(x -> x == JuMP.name(cref), constraint_names(owner_model(cref)))
    if length(idx) == 1
        delete!(object_dictionary(owner_model(cref)), Symbol(JuMP.name(cref)))
    end
    constraint_names(owner_model(cref))[index(cref)] = name
    return nothing
end

function JuMP.constraint_by_name(m::GMPModel, name::String)
    idx = findall(x -> x == name, constraint_names(m))
    if isempty(idx)
        return nothing
    elseif length(idx) > 1
        throw(ErrorException("Multiple constraints have the name $name."))
    else
        c = object_dictionary(m)[Symbol(name)]        
        return GMPConstraintRef(m, first(idx), JuMP.shape(c))
    end
end

function JuMP.show_constraints_summary(io::IO, m::GMPModel)
    n = length(gmp_constraints(m))
    print(io, "Constraint", JuMP._plural(n), ": ", n)
end

function JuMP.constraints_string(print_mode, m::GMPModel)
    strings = String[]
    # Sort by creation order, i.e. MOI.ConstraintIndex value
    constraints = sort(collect(gmp_constraints(m)), by = c -> first(c))
    for (index, constraint) in constraints
        push!(strings, JuMP.constraint_string(print_mode, constraint))
    end
    return strings
end

# Objective
function JuMP.set_objective(m::GMPModel, sense::MOI.OptimizationSense,
                            f::AbstractMomentExpressionLike)
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

function JuMP.show_objective_function_summary(mode, m::GMPModel)
    println(mode, "Objective function type: ",
            JuMP.objective_function_type(m))
end

function JuMP.objective_function_string(mode, m::GMPModel)
    return JuMP.function_string(mode, objective_function(m))
end

# Show
function JuMP.show_backend_summary(io::IO, model::GMPModel)
    println(io)
    println(io, "Approxmation mode: ", approximation_mode(model))
    println(io, "Maximum degree of data: ", model_degree(model))
    println(io, "Degree for approximation ", approximation_degree(model))
    # The last print shouldn't have a new line
    print(io, "Solver for approximation: ", solver_name(model))
end

# get value

function _is_approximated(vref::GMPVariableRef)
    if haskey(approx_vrefs(owner_model(vref)), index(vref))
        return true 
    else
        return false
    end
end

"""
    value(vref::GMPVariableRef)

Returns the approximation sequence of the referenced GMPObject.
"""
function JuMP.value(vref::GMPVariableRef)
    @assert _is_approximated(vref) "$vref has not been approximated yet"
    return approx_vrefs(owner_model(vref))[index(vref)].value 
end

export primal_justification
function primal_justification(vref::GMPVariableRef)
    @assert _is_approximated(vref) "$vref has not been approximated yet"

    return approx_vrefs(owner_model(vref))[index(vref)].p_just
end

export dual_justification
function dual_justification(vref::GMPVariableRef)
    @assert _is_approximated(vref) "$vref has not been approximated yet"    
    return approx_vrefs(owner_model(vref))[index(vref)].d_just
end

"""
approximation(vref::GMPVariableRef, m::AbstractPolynomialLike)

Return the approximation of the ⟨vref, m⟩. 
"""
function approximation(vref::GMPVariableRef, m::MP.AbstractPolynomialLike)
    @assert MP.variables(m) ⊆ MP.variables(vref) "$vref does not act on $m."
    @assert _is_approximated(vref) "$vref has not been approximated yet"
    as = approximation_sequence(owner_model(vref), index(vref)).dict
    if haskey(as, m)
        return value(as[m])
    else
        coefs, basis = MB.change_basis(m, approx_basis(vref_object(vref)))
        return sum( c*value(as[b]) for (c,b) in zip(coefs, basis)) 
    end
end

function approximation(vref::GMPVariableRef, m::Number)
    return approximation(vref,  m*_mono_one(variables(vref_object(vref))))
end

JuMP.value(me::MomentExpr) = sum(approximation(m, c) for (c, m) in me)

export residual
residual(cref::GMPConstraintRef) = approx_crefs(owner_model(cref))[index(cref)].p_res

function JuMP.dual(cref::GMPConstraintRef, ::Union{MomentConstraintShape, MeasureConstraintShape})
    return approx_crefs(owner_model(cref))[index(cref)].dual
end
