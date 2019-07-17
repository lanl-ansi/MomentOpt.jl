@testset "Macro Test" begin
    @testset "Model" begin
        gmp = GMPModel()
        @polyvar x y
        @test typeof(@measure gmp μ [x,y]) ==Vector{Measure}
        μ = gmp.measures[1]
        @measure(gmp, ν, [x,y], certificate=SOSCone())
        λ = Measure("λ",[x], support = @set(x^2-1>=0))
        @test typeof(@measure gmp λ) == Vector{Measure}
        @test typeof(@mconstraint gmp MomCons(Mom(μ,x^2)+π*Mom(y,ν),:eq,sqrt(2)))== Vector{MomCon{Polynomial{true,Float64},Float64}}
        @test typeof(@mconstraints gmp [MomCon(Mom(x^i*y, μ), :eq, i) for i=1:4])== Vector{MomCon{Polynomial{true,Float64},Float64}}
        @test typeof(@mconstraints gmp MomCons([Mom(x*i,λ)-Mom(y^i,μ) for i=1:10],:eq, 0)) == Vector{MomCon{Polynomial{true,Float64},Float64}}
        @test typeof(@objective gmp Max Mom(μ,x)+Mom(ν,y)) == MomExpr{Polynomial{true,Float64}}
    end
end
