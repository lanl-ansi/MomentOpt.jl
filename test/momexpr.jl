# test Mom, MomExpr
# test +,-,*,/

@testset "Moment Arith" begin
    @testset "Mom" begin
        @polyvar x y
        μ = Measure("μ", [x,y])
        @test Mom(1,μ) isa Mom{Int}
        @test LinearAlgebra.dot(x,μ) isa Mom{PolyVar{true}}
        @test Mom(μ,x*y) isa Mom{Monomial{true}}
        @test Mom(1,μ)+Mom(μ,x)	isa MomExpr{Polynomial{true,Int}}
        @test Mom(1.5*y,μ)+Mom(μ,x) isa MomExpr{Polynomial{true,Float64}}
        @test Mom(1,μ)-Mom(0.5,μ) isa MomExpr{Float64}
        @test Mom(1,μ)-Mom(μ,0.5*x) isa MomExpr{Polynomial{true,Float64}}
        @test 1.5*Mom(1,μ) isa Mom{Float64}
        @test Mom(1,μ)*3 isa Mom{Int}
        @test LinearAlgebra.dot(μ,1)/3 isa Mom{Float64}
        @test 1.5*Mom(y,μ) isa Mom{Term{true,Float64}}
        @test Mom(y,μ)*3 isa Mom{Term{true,Int}}
        @test Mom(y,μ)/3 isa Mom{Term{true,Float64}}
        @test [Mom(1,μ),LinearAlgebra.dot(μ,x)] isa Vector{Mom{Term{true,Int}}}
        @test [Mom(1,μ),LinearAlgebra.dot(1.5,μ)] isa Vector{Mom{Float64}}
        @test Mom.(μ,[1,x,x+y/3]) isa Vector{Mom{Polynomial{true,Float64}}}
        @polyvar z
        @test !MomentOpt.compatible(μ, z)
        ν = Measure("ν",[x])
        ERR = ErrorException("ν does not act on y")
        try Mom(ν,y)
        catch err
            @test err == ERR
        end
        ERR = ErrorException("The definition of a moment requires a measure and an integrand.")
        try Mom(y)
        catch err
            @test err == ERR
        end

        m1 = Mom(1,μ)
        @test MomentOpt.montype(m1) == Int
        m2 = Mom(μ, x*y)
        @test typeof([m1, m2]) == Vector{Mom{Term{true,Int}}}
    end

    @testset "MomExpr" begin
        @polyvar x y
        μ = Measure("μ", [x,y])
        ν = Measure("ν",[x])

        @test Mom.(μ, [1,x])-Mom.(μ,[y,2/3*y]) isa Vector{MomExpr{Polynomial{true,Float64}}}
        @test Mom.(μ, [1,x^2])+Mom.(μ,[y,y]) isa Vector{MomExpr{Polynomial{true,Int}}}
        @test MomExpr(1,μ) isa MomExpr{Int}
        @test MomExpr(μ,x) isa MomExpr{PolyVar{true}}
        @test MomExpr(OrderedDict{Measure,Int}(μ=>1)) isa MomExpr{Int}
        @test MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+π)) isa MomExpr{Polynomial{true,Float64}}
        @test MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>x^2+π)) isa MomExpr{Polynomial{true,Float64}}
        me1 = MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+π))
        me2 = MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>x^2+π))
        @test 3*me1 isa MomExpr{Polynomial{true,Float64}}
        @test me1+2*me2 isa MomExpr{Polynomial{true,Float64}}
        @test me1-me2 isa MomExpr{Polynomial{true,Float64}}
        @test MomExpr.([Mom(x+y,μ),Mom(1.5*x,μ)]) isa Vector{MomExpr{Polynomial{true,Float64}}}
        mev = MomExpr.([Mom(x+y,μ),Mom(1.5*x,μ)])
        @test mev+2*mev isa Vector{MomExpr{Polynomial{true,Float64}}}

        ERR = ErrorException("Measures and moments are not compatible.")
        try MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>y^2+π))
        catch err
            @test err==ERR
        end
        me3 = MomExpr(OrderedDict{Measure,Polynomial{true,Int}}(μ=>x^2+1))
        @test MomentOpt.montype(me3) == Polynomial{true,Int}
        @test promote_type(typeof(me1),typeof(me2)) == MomExpr{Polynomial{true,Float64}}
    end

    @testset "AffMomExpr" begin
        @polyvar x y
        μ = Measure("μ", [x,y])
        ν = Measure("ν",[x])
        m = Mom(y,μ)
        me = (Mom(1.5*y,μ)- Mom(1.5*x,ν))/3
        ae1 = AffMomExpr(me,1)
        @test m+1 isa AffMomExpr{PolyVar{true},Int}
        @test me+1 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1+1 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1+m isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1-me isa AffMomExpr{Polynomial{true,Float64},Int}
        ae2 = AffMomExpr(m,1)
        @test sprint(show, ae2) == "AffMomExpr{PolyVar{true},Int64}(⟨μ, y⟩, 1)" 
        @test ae1-ae2 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1+2*ae2 isa AffMomExpr{Polynomial{true,Float64},Int}

        @test promote_type(typeof(ae1),typeof(ae2)) == AffMomExpr{Polynomial{true,Float64},Int64}


    end
end
