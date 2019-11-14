@testset "Model Test" begin 
    m = GMPModel()
    @test Base.broadcastable(m) isa Base.RefValue{GMPModel}
    @test JuMP.object_dictionary(m) isa Dict
    @test JuMP.constraint_type(m) == ConstraintRef{GMPModel, Int64, MomentOpt.MomConShape}
    @polyvar x
	μ = Measure("μ", [x])
    @test add_measure!(m, μ) == Measure[μ]
    @test add_measure!(m, "ν", [x]) isa Vector{Measure}
    ν = measures(m)[2]
    @test measures(m) == [μ, ν]
    @test isempty(constraints(m))
    add_constraint(m, MomCon(Mom(μ,1), MOI.EqualTo(1)), "con")
    ϕ = Measure("ϕ", [x])
    @test set_objective(m, MOI.MIN_SENSE, Mom(μ, 1)) isa MomObj{Int}
    ERR = ErrorException("The model does not involve Measure{PolyVar{true},FullSpace,NonnegPolyInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle}}[ϕ]")
    try set_objective(m, MOI.MIN_SENSE, Mom(ϕ, 1))
       catch err
       @test err == ERR
    end 
end
