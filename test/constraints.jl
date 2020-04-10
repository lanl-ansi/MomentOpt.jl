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
    @polyvar x y
    m = GMPModel()
    @variable m μ Meas([x,y])
    @variable m ν Meas([x,y], basis = ChebyshevBasis)
    @variable m σ Meas([x])

    @test MOI.constant(MO.EqualToMeasure(ZeroMeasure())) isa ZeroMeasure
    @test JuMP.reshape_vector(μ + ν , MeasureConstraintShape()) == μ + ν
    @test  JuMP.reshape_set(EqualToMeasure(ZeroMeasure), MeasureConstraintShape()) isa ZeroMeasure()
    @test MO.same_basis(MO.GMPVariable[]) isa Nothing 
    @test MO.same_basis([μ, σ]) isa Nothing
    @test MO.same_basis([μ, ν ]) isa Nothing

    mc = MO.MeasureConstraint(μ + ν, MO.EqualToMeasure(ZeroMeasure())) 
    @test sprint(show, mc) == "μ + ν = 0"
    @test JuMP.jump_function(mc) == μ + ν
    @test JuMP.moi_set(mc) == MO.EqualTo(ZeroMeasure())
    @test JuMP.shape(mc) isa MeasureConstraintShape

    _err = x -> ""
    @test JuMP.build_constraint(_err, μ + ν - 0, MOI.EqualTo(0)) isa MO.MeasureConstraint
    @test_throws AssertionError JuMP.build_constraint(_err, μ + ν - 0, MOI.EqualTo(1))
    @test_throws AssertionError JuMP.build_constraint(_err, μ + ν - 1, MOI.EqualTo(0))

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
