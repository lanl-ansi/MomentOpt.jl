@testset "measexpr" begin
    @polyvar x y
    m = GMPModel()
    @measure m μ[1:2] [x]
    @measure m ν [x,y]
    me1 = MeasureExpr([0.5, 0.5], μ)
    
    @test sprint(show, me1) == "0.5μ[1] + 0.5μ[2]"
    @test sprint(show, -μ[1]) == " - μ[1]"
    @test coefficients(me1) == [0.5, 0.5]   
    @test MomentOpt.measures(me1) == μ 
    @test poly_variables(me1) == [x]

    me2 = MeasureExpr([0.5, -0.5], μ)
    @test sprint(show, me2) == "0.5μ[1] - 0.5μ[2]"

    @test_throws AssertionError MeasureExpr(ones(3)/3, μ)
    @test_throws ArgumentError MeasureExpr(ones(3)/3, [μ...,ν])

    @measure m ρ [x]
    @test [me1, ρ] isa Vector{MeasureExpr{T}} where T <: AbstractFloat    
    @test [me1, MeasureExpr([1], [ρ])] isa Vector{MeasureExpr{T}} where T <: AbstractFloat

    @test coefficients(-me1) == [-0.5,-0.5] 
    @test MomentOpt.measures(-me1) == MomentOpt.measures(me1)
    @test coefficients(me1*2) == [1.0, 1.0]
    @test coefficients(me1/2) == [0.25, 0.25]

    @test length(coefficients(me1 + me2)) == 1
    @test length(coefficients(me1 - me2)) == 1
    
    @test sum(μ) isa MeasureExpr{Int}
    @test sum(μ) - μ[1] == μ[2]

end

