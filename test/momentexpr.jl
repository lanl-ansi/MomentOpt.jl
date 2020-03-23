# test MomentExpr
# test +,-,*,/

    @testset "MomExpr" begin
        @polyvar x y
        μ = Measure("μ", [x,y])
        ν = Measure("ν",[x])

        @test Mom.(μ, [1,x])-Mom.(μ,[y,using) isa Vector{MomExpr{Polynomial{true,Float64}}}
        @test Mom.(μ, [1,x^2])+Mom.(μ,[y,y]) isa Vector{MomExpr{Polynomial{true,Int}}}
        @test MomExpr(1,μ) isa MomExpr{Int}
        @test MomExpr(μ,x) isa MomExpr{PolyVar{true}}
        @test MomExpr(OrderedDict{Measure,Int}(μ=>1)) isa MomExpr{Int}
        @test MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+1.2)) isa MomExpr{Polynomial{true,Float64}}
        @test MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>x^2+1.2)) isa MomExpr{Polynomial{true,Float64}}
        me1 = MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+1.2))
        me2 = MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>x^2+1.2))
        @test 3*me1 isa MomExpr{Polynomial{true,Float64}}
        @test me1+2*me2 isa MomExpr{Polynomial{true,Float64}}
        @test me1-me2 isa MomExpr{Polynomial{true,Float64}}
        @test MomExpr.([Mom(x+y,μ),Mom(1.5*x,μ)]) isa Vector{MomExpr{Polynomial{true,Float64}}}
        mev = MomExpr.([Mom(x+y,μ),Mom(1.5*x,μ)])
        @test mev+2*mev isa Vector{MomExpr{Polynomial{true,Float64}}}

        ERR = ErrorException("Measures and moments are not compatible.")
        try MomExpr(OrderedDict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>y^2+1.2))
        catch err
            @test err==ERR
        end
        me3 = MomExpr(OrderedDict{Measure,Polynomial{true,Int}}(μ=>x^2+1))
        @test MomentOpt.montype(me3) == Polynomial{true,Int}
        @test Base.promote_rule(typeof(me1), typeof(Mom(1,μ))) == MomExpr{Polynomial{true,Float64}}
        @test promote_type(typeof(me1),typeof(me2)) == MomExpr{Polynomial{true,Float64}}
    end

    @testset "AffMomExpr" begin
        @polyvar x y
        μ = Measure("μ", [x,y])
        ν = Measure("ν",[x])
        m = Mom(y,μ)
        me = (Mom(1.5*y,μ)- Mom(1.5*x,ν))/3

        @test MomentOpt.add_mom_type([Mom(1,ν), Mom(μ,2)]) == Int
        @test MomentOpt.add_mom_type([m,me]) == Polynomial{true,Float64}
        @test me - 1 isa AffMomExpr{Polynomial{true,Float64},Int64}
        @test 1 - me isa AffMomExpr{Polynomial{true,Float64},Int64} 
        
        ae1 = AffMomExpr(me,1)
        
        @test m+1 isa AffMomExpr{PolyVar{true},Int}
        @test ae1 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1+1 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1+m isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1-me isa AffMomExpr{Polynomial{true,Float64},Int}        
        @test 1+m isa AffMomExpr{PolyVar{true},Int}
        @test 1+ae1 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test m+ae1 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test me-ae1 isa AffMomExpr{Polynomial{true,Float64},Int}
        

        ae2 = AffMomExpr(m,1)
        @test sprint(show, ae2) == "⟨μ, y⟩ + 1"
        @test sprint(show, -ae2) == "⟨μ, -y⟩ - 1"
        @test sprint(show, AffMomExpr(m,0)) == "⟨μ, y⟩"
        @test ae1-ae2 isa AffMomExpr{Polynomial{true,Float64},Int}
        @test ae1+2*ae2 isa AffMomExpr{Polynomial{true,Float64},Int}

        @test promote_rule(AffMomExpr{Polynomial{true,Int}, Int}, AffMomExpr{Polynomial{true,Float64},Float64})== AffMomExpr{Polynomial{true,Float64},Float64}
        @test promote_rule(typeof(ae1),typeof(m)) == AffMomExpr{Polynomial{true,Float64},Int64}
        @test promote_rule(typeof(ae1),typeof(me)) == AffMomExpr{Polynomial{true,Float64},Int64}

        
    end
end
