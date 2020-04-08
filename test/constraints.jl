@testset "MomentConstraint" begin
# JuMP.reshape_vector(expr::MomentExpr, ::MomentConstraintShape) 
# JuMP.reshape_set(set::MOI.AbstractScalarSet, ::MomentConstraintShape)
# MomentConstraint(ae::AffMomentExpr, set)
# JuMP.jump_function(con::MomentConstraint)
# JuMP.moi_set(con::MomentConstraint)
# JuMP.shape(con::MomentConstraint)
# JuMP.function_string(mode, mc::MomentConstraint)
# Base.show(io::IO, con::MomentConstraint)
# JuMP.build_constraint(_error::Function, ae::AffMomentExpr, set::MOI.AbstractScalarSet)
end

@testset "MeasureConstraint" begin
# MOI.constant(set::EqualToMeasure) 
# JuMP.in_set_string(print_mode, set::EqualToMeasure)
# JuMP.reshape_vector(expr::ObjectExpr, ::MeasureConstraintShape)
# JuMP.reshape_set(set::EqualToMeasure, ::MeasureConstraintShape)
# same_basis(mv::Vector{GMPVariableRef})
#
# JuMP.jump_function(con::MeasureConstraint)
# JuMP.moi_set(con::MeasureConstraint)
# JuMP.shape(con::MeasureConstraint)
# JuMP.function_string(mode, mc::MeasureConstraint)
# Base.show(io::IO, con::MeasureConstraint)
# JuMP.build_constraint(_error::Function, ae::AffObjectExpr, s::MOI.EqualTo)
end

@testset "GMPConstraintRef" begin
# Base.iszero(::GMPConstraintRef) 
# JuMP.isequal_canonical(v::GMPConstraintRef, w::GMPConstraintRef)
# Base.:(==)(v::GMPConstraintRef, w::GMPConstraintRef)
# JuMP.owner_model(con_ref::GMPConstraintRef)
# JuMP.shape(cref::GMPConstraintRef)
# JuMP.index(cref::GMPConstraintRef)
end
