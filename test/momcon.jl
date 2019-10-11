@testset "MomCon Test" begin

    @testset "MomObj" begin
        @test Base.broadcastable(MomentOpt.MomConShape()) isa Base.RefValue{MomentOpt.MomConShape}(MomentOpt.MomConShape())
        @polyvar x y
        μ = Measure("μ", [x,y])
        ν = Measure("ν",[x])
        @test MomObj(MOI.MAX_SENSE, μ, 1) isa MomObj{Int}
        @test MomObj(MOI.MIN_SENSE, x^2, ν) isa MomObj{Monomial{true}}
        obj = MomObj(MOI.MAX_SENSE, Mom(μ, y) + Mom(ν, x))
        @test typeof(obj) == MomObj{Polynomial{true, Int}}
        @test ν in measures(obj)
        @test μ in measures(obj)
    end

    @testset "MomCon" begin
        @polyvar x[1:3]
        μ = Measure("μ", x)
        ν = Measure("ν", x)
        @test typeof(MomCon(Mom(μ, x[1]), MOI.EqualTo(1))) == MomCon{PolyVar{true}}
        @test typeof(MomCon(Mom(μ, x[1]) + Mom(ν,x[2]), MOI.LessThan(0)))== MomCon{Polynomial{true, Int}}
        @test typeof(MomCon(Mom(μ, 1.5), MOI.EqualTo(0))) == MomCon{Float64}
        @test typeof([MomCon(Mom(x[i],μ) + Mom(ν,x[i]), MOI.EqualTo(0)) for i =1:3])== Vector{MomCon{Polynomial{true, Int}}}
        @test typeof(MomCon.([Mom(x[i],μ) + Mom(ν,x[i]) for i = 1:3], MOI.EqualTo(0)))== Vector{MomCon{Polynomial{true, Int}}}

        mc = MomCon(Mom(μ, x[1]), MOI.EqualTo(1))
        mc1 = MomCon(Mom(μ, 1.0*x[1]), MOI.EqualTo(1))
        
        @test sprint(show, mc) == "⟨μ, x[1]⟩ = 1" 
        @test MomentOpt.constant(mc) == 1
        @test JuMP.jump_function(mc) isa MomExpr
        @test JuMP.shape(mc) == MomentOpt.MomConShape()
        @test JuMP.reshape_vector(mc.func, MomentOpt.MomConShape()) == mc.func
        @test JuMP.reshape_set(mc.set, MomentOpt.MomConShape()) == mc.set

        @test Base.promote_rule(typeof(mc),typeof(mc1)) == MomCon{Term{true,Float64}}
        mm = [mc, mc1]
        @test mm isa Array{MomCon{Term{true, Float64}},1}
        @test measures(mm) == Measure[μ]

        @test MomentOpt.constant(MOI.EqualTo(1)) == 1
        @test MomentOpt.constant(MOI.LessThan(1)) == 1
        @test MomentOpt.constant(MOI.GreaterThan(1)) == 1

    end
end
