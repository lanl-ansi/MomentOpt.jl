@testset "MomCon Test" begin
    @testset "MomObj" begin
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
        @test typeof(MomCon(Mom(μ,x[1]) + Mom(ν,x[2]), MOI.LessThan(0)))== MomCon{Polynomial{true, Int}}
        @test typeof(MomCon(Mom(μ, 1.5), MOI.EqualTo(0))) == MomCon{Float64}
        @test typeof([MomCon(Mom(x[i],μ) + Mom(ν,x[i]), MOI.EqualTo(0)) for i =1:3])== Vector{MomCon{Polynomial{true, Int}}}
        @test typeof(MomCon.([Mom(x[i],μ) + Mom(ν,x[i]) for i = 1:3], MOI.EqualTo(0)))== Vector{MomCon{Polynomial{true, Int}}}        
    end
end
