@testset "Constraints" begin
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

        @test MOI.constant(MO.EqualToMeasure(ZeroMeasure([x, y]))) isa AnalyticMeasure
        @test JuMP.reshape_vector(μ + ν , MO.MeasureConstraintShape()) == μ + ν
        @test  JuMP.reshape_set(MO.EqualToMeasure(ZeroMeasure([x])), MO.MeasureConstraintShape()) isa MO.EqualToMeasure

        mc = MO.MeasureConstraint(μ + ν, MO.EqualToMeasure(ZeroMeasure([x,y]))) 
        @test sprint(show, mc) == "μ + ν = AnalyticMeasure"
        @test JuMP.jump_function(mc) == μ + ν
        @test JuMP.moi_set(mc) isa MO.EqualToMeasure
        @test JuMP.shape(mc) isa MO.MeasureConstraintShape

        _err = x -> ""
        @test JuMP.build_constraint(_err, μ + ν, MO.EqualToMeasure(ZeroMeasure([x, y]))) isa MO.MeasureConstraint
     end

    @testset "GMPConstraintRef" begin
        @polyvar x y
        m = GMPModel()
        @variable m μ Meas([x,y])
        @constraint m c1 Mom(1, μ) == 1
        @constraint m c2 μ == ZeroMeasure([x,y]) 
        @test !iszero(c1)
        @test !JuMP.isequal_canonical(c1, c2)
        @test JuMP.isequal_canonical(c1, c1) 
        @test JuMP.isequal_canonical(c2, c2) 
        @test JuMP.owner_model(c1) == m
        @test JuMP.shape(c1) isa MO.MomentConstraintShape
        @test JuMP.shape(c2) isa MO.MeasureConstraintShape
        @test JuMP.index.([c1, c2]) == [1, 2] 
    end

    @testset "Add/DeleteConstraint" begin
        @polyvar x y
        m = GMPModel()
        @variable m μ Meas([x,y])
        @variable m ν Meas([x,y], basis = ChebyshevBasis)
        @constraint m c1 Mom(1, μ) == 1

        ζ = ZeroMeasure([x,y]) 
        @constraint m c2 μ == ζ
        @test JuMP.is_valid(m, c1)
        @test sprint(show, c1) == "c1 : ⟨1, μ⟩ = 1.0"
        @test sprint(show, c2) == "c2 : μ = AnalyticMeasure"

        @test sprint(show, JuMP.constraint_object(c2)) == "μ = AnalyticMeasure"

        @test JuMP.jump_function(c2) == MO.ObjectExpr(1, μ)
        @test JuMP.moi_set(c2) isa MO.EqualToMeasure
        @test JuMP.name(c1) == "c1"
        @test JuMP.name(c2) == "c2"

        JuMP.set_name(c1, "c2")
        @test JuMP.name(c1) == "c2"

        @test JuMP.constraint_by_name(m, "c1") isa Nothing
        @test_throws ErrorException JuMP.constraint_by_name(m, "c2")

        JuMP.delete(m, c1)
        @test JuMP.constraint_by_name(m, "c2") == c2

        @test JuMP.num_constraints(m, MO.MeasureConstraint) == 1
        @test JuMP.num_constraints(m, MO.MomentConstraint) == 0

    end

end
