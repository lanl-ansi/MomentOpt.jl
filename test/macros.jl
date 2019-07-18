@testset "Macro Test" begin
    @testset "Model" begin
        gmp = GMPModel()
        @polyvar x y
        @test typeof(@measure gmp μ [x, y]) == Vector{Measure}
        μ = gmp.measures[1]
        @measure(gmp, ν, [x,y], certificate = SOSCone())
        λ = Measure("λ", [x], support = @set(x^2-1>=0))
        @test typeof(@measure gmp λ) ==  Vector{Measure}
        @test typeof(@objective gmp Max Mom(μ,x) + Mom(ν,y)) == MomExpr{Polynomial{true,Float64}}
        @test typeof(@constraint gmp Mom(μ,y) == 1) == ConstraintRef{GMPModel,Int64,MomentOpt.MomConShape}
        @test typeof(@constraint gmp c[i=1:5] Mom(λ,x^(i-1)) == 1/i) == Array{ConstraintRef{GMPModel,Int64,MomentOpt.MomConShape},1}
    end
end



