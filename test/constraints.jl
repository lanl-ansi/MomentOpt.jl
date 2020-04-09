@testset "MomentConstraint" begin
    @polyvar x y
    m = GMPModel()
    @variable m μ Meas([x,y])
    @variable m ν Meas([x,y], basis = ChebyshevBasis)
    @variable m σ Meas([x])

    oe = μ + ν
    me = Mom(x*y + x, oe)
    @test JuMP.reshape_vector(me, MO.MomentConstraintShape()) == me
    @test JuMP.reshape_set(MOI.EqualTo(0), MO.MomentConstraintShape()) == MOI.EqualTo(0)

    con1 = MO.MomentConstraint(me, MOI.EqualTo(0))
    con2 = MO.MomentConstraint(me+1, MOI.EqualTo(1))

    @test JuMP.jump_function(con1) == JuMP.jump_function(con2) 
    @test JuMP.moi_set(con1) == JuMP.moi_set(con2) 
    @test JuMP.shape(con1) isa MO.MomentConstraintShape
    @test sprint(show, con1) == sprint(show, con2)

    _err = x -> ""
    con3 = JuMP.build_constraint(_err, me+1, MOI.EqualTo(1))
    @test JuMP.jump_function(con3) == JuMP.jump_function(con2) 
    @test JuMP.moi_set(con3) == JuMP.moi_set(con2) 

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

@testset "Add/DeleteConstraint" begin

end
