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
    objective_function::AbstractGMPExpressionLike
    max_degree::Int
    function ModelData()
        new(Dict{Int, AbstractGMPVariable}(),
            Dict{Int, AbstractGMPConstraint}(),
            MOI.FEASIBILITY_SENSE,
            GMPEmptyExpr(), 0)
    end
end

# define types informing about the mode and status of the model.
"""
    AbstractAproximationMode

Defines whether the primal or the dual of a problem should be used for a conic approximation.
"""

export NO_APPROXIMATION_MODE, PRIMAL_RELAXATION_MODE, DUAL_RELAXATION_MODE, PRIMAL_STRENGTHEN_MODE, DUAL_STRENGTHEN_MODE

abstract type AbstractApproximationMode end
struct NO_APPROXIMATION_MODE <: AbstractApproximationMode end 
struct PRIMAL_STRENGTHEN_MODE <: AbstractApproximationMode end 
struct PRIMAL_RELAXATION_MODE <: AbstractApproximationMode end 
struct DUAL_STRENGTHEN_MODE <: AbstractApproximationMode end
struct DUAL_RELAXATION_MODE <: AbstractApproximationMode end
mutable struct ApproximationInfo
    mode::AbstractApproximationMode
    degree::Int
    solver::Union{Nothing, Function}
    approx_vrefs::Dict{Int, Any}
    approx_crefs::Dict{Int, Any}
    function ApproximationInfo()
        new(PRIMAL_RELAXATION_MODE(), 0, nothing, Dict{GMPVariableRef, Any}(), Dict{GMPConstraintRef, Any}())
    end
end

# define GMPmodel <: JuMP.AbstractModel
# The implementation follows basically the example in JuMP/test/JuMPExtension.jl
export GMPModel
"""
    GMPModel

This JuMP.AbstractModel is able to model linear obtimization problems where the variables are abstract measures or continuous functions.
"""
mutable struct GMPModel <: JuMP.AbstractModel
    model_info::ModelInfo
    model_data::ModelData    
    approximation_info::ApproximationInfo
    approximation_model::Union{Nothing, JuMP.Model}
    function GMPModel()
        new(ModelInfo(), ModelData(), ApproximationInfo(), nothing)
    end
end
Base.broadcastable(model::GMPModel) = Ref(model)
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
    return termination_status(approximation_model(model))
end

function JuMP.set_optimizer(model::GMPModel, optimizer::Function)
    approximation_info(model).solver = optimizer
    return nothing
end


function update_degree(m::GMPModel, degree::Int)
    if model_data(m).max_degree < degree
        model_data(m).max_degree = degree
        approximation_info(m).degree = maximum([approximation_info(m).degree, degree])
    end
end

export set_approximation_degree
"""
    set_approximation_degree(model::GMPModel, degree::Int)

Manually set the degree of the approximation to a custom value. The value cannot be less than the degree of the data. 
"""
function set_approximation_degree(model::GMPModel, degree::Int)
    if approximation_info(model).degree > degree && model_data(model).max_degree > degree
        approximation_info(model).degree = model_data(model).max_degree
        @warn "Requested approximation degree $degree is too low to cover all data. The approximation degree has been set to the minimal value possible and now is $(model_data(model).max_degree)  "
    else
        approximation_info(model).degree = degree
    end
    return nothing
end

export set_approximation_mode
"""
    set_approximation_mode(m::GMPModel, mode::AbstractApproximationMode) 

Set the approximation mode of m to mode. Possible values for mode are:
  * NO_APPROXIMATION_MODE()
  * PRIMAL_STRENGTHEN_MODE() 
  * PRIMAL_RELAXATION_MODE() (default)
  * DUAL_STRENGTHEN_MODE()
  * DUAL_RELAXATION_MODE() 
  """
function set_approximation_mode(m::GMPModel, mode::AbstractApproximationMode) 
    approximation_info(m).mode = mode
    return m
end

# referencing variables
JuMP.variable_type(::GMPModel) = GMPVariableRef

function GMPVariableRef(m::GMPModel, v::AbstractGMPVariable)
    model_info(m).vct += 1
    return GMPVariableRef(m, model_info(m).vct, variable_type(v))
end

JuMP.name(vref::GMPVariableRef) = variable_names(vref.model)[vref.index]

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
        v = object_dictionary(m)[Symbol(name)]
        return GMPVariableRef(m, first(idx), typeof(v))
    end
end

function JuMP.is_valid(m::GMPModel, vref::GMPVariableRef)
    return (m === owner_model(vref) &&  index(vref) in keys(gmp_variables(m)))
end

function JuMP.add_variable(m::GMPModel, v::AbstractGMPVariable, name::String = "")
    vref = GMPVariableRef(m, v)
    gmp_variables(m)[index(vref)] = v
    push!(variable_names(m), name)
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
    return MOI.get(vref_object(vref), Variables())
end

export support
"""
    support(vref::GMPVariableRef)

Retruns the support a measure is defined on. 
"""
function support(vref::GMPVariableRef)
    return MOI.get(vref_object(vref), BSASet())
end

export domain
"""
    domain(vref::GMPVariableRef)

Retruns the domain a continuous is defined on. 
"""
function domain(vref::GMPVariableRef)
    return MOI.get(vref_object(vref), BSASet())
end

export approx_basis
"""
    approx_basis(vref::GMPVariableRef)

Returns the basis used to approximate vref.
"""
function approx_basis(vref::GMPVariableRef)
    return MOI.get(vref_object(vref), ApproximationBasis())
end

# referencing constraints 
JuMP.constraint_type(::GMPModel) = GMPConstraintRef

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

function JuMP.add_constraint(m::GMPModel, c::ScalarConstraint{MomentOpt.ObjectExpr{S}, MOI.EqualTo{T}}, name::String = "") where {S, T}
    @assert(MOI.constant(moi_set(c)) == 0, "Affine measure constraints are not supported.")    
    con = MeasureConstraint(jump_function(c), EqualToMeasure(ZeroMeasure()))
    return add_constraint(m, con, name)
end


function JuMP.delete(m::GMPModel, cref::GMPConstraintRef)
    @assert is_valid(m, cref)
    delete!(gmp_constraints(m), index(cref))
    delete!(object_dictionary(m), index(cref))
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

function JuMP.show_constraints_summary(io::IO, m::GMPModel)
    n = length(gmp_constraints(m))
    print(io, "Constraint", JuMP._plural(n), ": ", n)
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

# Objective
function JuMP.set_objective(m::GMPModel, sense::MOI.OptimizationSense,
                            f::AbstractGMPExpressionLike)
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

function JuMP.show_objective_function_summary(io::IO, m::GMPModel)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(m))
end

function JuMP.objective_function_string(print_mode, m::GMPModel)
    return JuMP.function_string(print_mode, objective_function(m))
end


# Show
function JuMP.show_backend_summary(io::IO, model::GMPModel)
    println(io)
    println(io, "Approxmation mode: ", approximation_info(model).mode)
    println(io, "Maximum degree of data: ", model_data(model).max_degree)
    println(io, "Degree for approximation ", approximation_info(model).degree)
    # The last print shouldn't have a new line
    if approximation_info(model).solver isa Nothing
        print(io, "Solver for approximation: ")
    else
        print(io, "Solver for approximation: ", approximation_info(model).solver)
    end
end

# get value
export approximation
"""
    approximation(vref::GMPVariableRef)

Returns the approximation sequence of the referenced GMPObject.
"""
function approximation(vref::GMPVariableRef)
    return Dict(k => value(v) for (k, v) in approximation_sequence(owner_model(vref), index(vref)).dict)
end

"""
    approximation(vref::GMPVariableRef, m::AbstractPolynomialLike)

Return the approximation of the ⟨vref, m⟩. 
"""
function approximation(vref::GMPVariableRef, m::MP.AbstractPolynomialLike)
    #TODO check whether approximation is available
    as = approximation_sequence(owner_model(vref), index(vref)).dict
    if haskey(as, m)
        return value(as[m])
    else
        coefs, basis = MB.change_basis(m, MOI.get(vref_object(vref), ApproximationBasis()))
        return sum( c*value(as[b]) for (c,b) in zip(coefs, basis)) 
    end
end

function approximation(vref::GMPVariableRef, m::Number)
    return approximation(vref,  m*_mono_one(MOI.get(vref_object(vref), Variables())))
end


JuMP.value(me::MomentExpr) = sum(c*approximation(m, cv) for (cv, mv) in me for (c, m) in mv)
