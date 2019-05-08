@testset "MomCon Tests" begin
    @testset "MomObj" begin
        @polyvar x y
        μ = Measure("μ", [x,y])
        ν = Measure("ν",[x])
        @test typeof(MomObj(:Max, μ, 1)) == CMomObj{Int}
        @test typeof(MomObj(:Min, x^2,ν)) == MomObj{Monomial{true}}
        obj = MomObj(:max, Mom(μ, y)+Mom(ν, x))
        @test typeof(obj) == MomObj{Polynomial{true,Int}}
        @test ν in measures(obj)
	@test μ in measures(obj)
    end
    @testset "MomCon" begin
        @polyvar x[1:3]
        μ = Measure("μ", x)
        ν = Measure("ν", x)
        @test typeof(MomCon(Mom(μ, x[1]), :eq, 1)) == MomCon{PolyVar{true},Int}
        @test typeof(MomCon(Mom(μ,x[1])+Mom(ν,x[2]),:leq, 0))==MomCon{Polynomial{true,Int},Int}
        @test typeof(MomCon(Mom(μ, 1.5),:eq, 0)) == CMomCon{Float64,Int}
        @test typeof([MomCon(Mom(x[i],μ)+Mom(ν,x[i]),:eq,0) for i =1:3])== Vector{MomCon{Polynomial{true,Int},Int}}
    end
    @testset "MomCons" begin
        @polyvar x[1:3]
        μ = Measure("μ", x)
        ν = Measure("ν", x)
        @test typeof(MomCons(Mom(μ,1),:eq,0)) == CMomCon{Int,Int}
        momcon = MomCons([Mom(x[i]^i,μ) for i = 1:3],:eq, collect(1:3))
        @test typeof(momcon) == Vector{MomCon{Monomial{true},Int}}
        @test μ in measures(momcon)
	@test !(ν in measures(momcon))
        @test typeof(MomCons([Mom(μ,x'*x),Mom(ν,x'*x)],[:geq,:leq],[0,1])) == Vector{MomCon{Polynomial{true,Int},Int}}
    end

end
