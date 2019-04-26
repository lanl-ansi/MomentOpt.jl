# test  CMom, Mom, CMomExpr, MomExpr
# test +,-,*,/

@testset "Moment Arithmetic" begin
	@testset "Mom" begin
		@polyvar x y
		μ = Measure("μ", [x,y])
		@test Mom(1,μ) isa CMom{Int}
		@test Mom(μ,x) isa Mom{PolyVar{true}}
		@test Mom(μ,x*y) isa Mom{Monomial{true}}
		@test Mom(1,μ)+Mom(μ,x)	isa Mom{Polynomial{true,Int}}
		@test Mom(1.5*y,μ)+Mom(μ,x) isa Mom{Polynomial{true,Float64}}
		@test Mom(1,μ)-Mom(0.5,μ) isa CMom{Float64}
		@test Mom(1,μ)-Mom(μ,0.5*x) isa Mom{Polynomial{true,Float64}}
		@test 1.5*Mom(1,μ) isa CMom{Float64}
		@test Mom(1,μ)*3 isa CMom{Int}
		@test Mom(1,μ)/3 isa CMom{Float64}
		@test 1.5*Mom(y,μ) isa Mom{Term{true,Float64}}
		@test Mom(y,μ)*3 isa Mom{Term{true,Int}}
		@test Mom(y,μ)/3 isa Mom{Term{true,Float64}}
		@test [Mom(1,μ),Mom(μ,x)] isa Vector{Mom{Term{true,Int}}}
		@test [Mom(1,μ),Mom(1.5,μ)] isa Vector{CMom{Float64}}
		@test Mom(μ,[1,x,x+y/3]) isa Vector{Mom{Polynomial{true,Float64}}}
		ν = Measure("ν",[x])
		@test Mom(y,ν) == nothing
	end
	
	@testset "MomExpr" begin
		@polyvar x y
		μ = Measure("μ", [x,y])
		ν = Measure("ν",[x])
		@test Mom(μ, [1,x])-Mom(μ,[y,2/3*y]) isa Vector{MomExpr{Polynomial{true,Float64}}}
		@test Mom(μ, [1,x^2])+Mom(μ,[y,y]) isa Vector{MomExpr{Polynomial{true,Int}}}
		@test MomExpr(1,μ) isa CMomExpr{Int}
		@test MomExpr(μ,x) isa MomExpr{PolyVar{true}}
		@test MomExpr(Dict{Measure,Int}(μ=>1)) isa CMomExpr{Int}
		@test MomExpr(Dict{Measure,Polynomial{true,Float64}}(μ=>x^2+π)) isa MomExpr{Polynomial{true,Float64}}
		@test MomExpr(Dict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>x^2+π)) isa MomExpr{Polynomial{true,Float64}}
		me1 = MomExpr(Dict{Measure,Polynomial{true,Float64}}(μ=>x^2+π)) 
		me2 = MomExpr(Dict{Measure,Polynomial{true,Float64}}(μ=>x^2+y,ν=>x^2+π)) 
		@test 3*me1 isa MomExpr{Polynomial{true,Float64}}
		@test me1+2*me2 isa MomExpr{Polynomial{true,Float64}}
		@test me1-me2 isa MomExpr{Polynomial{true,Float64}}
		@test MomExpr([Mom(x+y,μ),Mom(1.5*x,μ)]) isa Vector{MomExpr{Polynomial{true,Float64}}}
		mev = MomExpr([Mom(x+y,μ),Mom(1.5*x,μ)])
		@test mev+2*mev isa Vector{MomExpr{Polynomial{true,Float64}}}
	end

	
end




