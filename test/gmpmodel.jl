@testset "GMPModel" begin
    @testset "Model" begin
# MO.model_info(model::GMPModel) = model.model_info
# MO.model_data(model::GMPModel) = model.model_data
# MO.approximation_info(model::GMPModel) = model.approximation_info
# MO.approximation_model(model::GMPModel) = model.approximation_model
# MO.gmp_variables(model::GMPModel) = model_data(model).variables
# MO.gmp_constraints(model::GMPModel) = model_data(model).constraints
# MO.variable_names(model::GMPModel) = model_info(model).variable_names
# MO.constraint_names(model::GMPModel) = model_info(model).constraint_names
# JuMP.object_dictionary(model::GMPModel) = model.model_info.obj_dict
# JuMP.num_variables(model::GMPModel) = length(gmp_variables(model))
# JuMP.termination_status(model::GMPModel)
# JuMP.variable_type(::GMPModel)
# 
# JuMP.set_optimizer(model::GMPModel, optimizer::Function)
# MO.update_degree(m::GMPModel, degree::Int)
# set_approximation_degree(model::GMPModel, degree::Int)
# set_approximation_mode(m::GMPModel, mode::AbstractApproximationMode) 
end
@testset "Variables" begin
# JuMP.name(vref::GMPVariableRef)

# JuMP.set_name(v::GMPVariableRef, s::String)
# JuMP.is_valid(m::GMPModel, vref::GMPVariableRef)
# GMPVariableRef(m::GMPModel, v::AbstractGMPVariable)
# JuMP.add_variable(m::GMPModel, v::AbstractGMPVariable, name::String = "")
# JuMP.delete(m::GMPModel, vref::GMPVariableRef)
# JuMP.variable_by_name(m::GMPModel, s::String)
# MO.vref_object(vref::GMPVariableRef)
# MP.variables(vref::GMPVariableRef)
#
# support(vref::GMPVariableRef)
# domain(vref::GMPVariableRef)
# approx_basis(vref::GMPVariableRef)
end
@testset "Constraints" begin
# JuMP.constraint_type(::GMPModel)
# JuMP.is_valid(m::GMPModel, cref::GMPConstraintRef)
# JuMP.constraint_string(print_mode, ref::GMPConstraintRef; in_math_mode = false)
# Base.show(io::IO, ref::GMPConstraintRef)
# JuMP.constraint_object(cref::GMPConstraintRef)
# JuMP.jump_function(cref::GMPConstraintRef)
# JuMP.moi_set(cref::GMPConstraintRef)
# JuMP.add_constraint(m::GMPModel, c::AbstractGMPConstraint, name::String = "")
# JuMP.add_constraint(m::GMPModel, c::ScalarConstraint{MomentOpt.ObjectExpr{S}, MOI.EqualTo{T}}, name::String = "") where {S, T}
# JuMP.delete(m::GMPModel, cref::GMPConstraintRef)
# JuMP.num_constraints(m::GMPModel, F::Type{<:JuMP.AbstractJuMPScalar}, S::Type{<:MOI.AbstractSet})
# MP.num_constraints(m::GMPModel, ::Type{<:Vector{F}}, S::Type{<:MOI.AbstractSet}) where F<:JuMP.AbstractJuMPScalar
# JuMP.name(cref::GMPConstraintRef) 
# JuMP.set_name(cref::GMPConstraintRef, name::String)
# JuMP.constraint_by_name(m::GMPModel, name::String)
# JuMP.show_constraints_summary(io::IO, m::GMPModel)
# JuMP.constraints_string(print_mode, m::GMPModel)

end
@testset "Objective" begin
# JuMP.set_objective(m::GMPModel, sense::MOI.OptimizationSense, f::AbstractGMPExpressionLike)
# JuMP.objective_sense(m::GMPModel) 
# JuMP.set_objective_sense(m::GMPModel, sense::MOI.OptimizationSense)
# JuMP.objective_function(m::GMPModel) 
# JuMP.objective_function_type(m::GMPModel) 
# JuMP.objective_function(m::GMPModel, FT::Type)
# JuMP.objective_function_string(print_mode, m::GMPModel)
end
@testset "Summary" begin
# JuMP.show_objective_function_summary(io::IO, m::GMPModel)
# JuMP.show_backend_summary(io::IO, model::GMPModel)
end

@testset "Get Value" begin
    # approximation(vref::GMPVariableRef)
    # approximation(vref::GMPVariableRef, m::MP.AbstractPolynomialLike)
    # approximation(vref::GMPVariableRef, m::Number)
    # JuMP.value(me::MomentExpr)
end

end

