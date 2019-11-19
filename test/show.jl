@testset "Show Test" begin
    gmp = GMPModel()
    @test sprint(show, gmp) == """
    GMPModel:
    Feasibility problem:
    Unknowns:
    OPTIMIZE_NOT_CALLED
    """
    @polyvar x y
    @measure gmp μ [x,y]
    @measure gmp ν [x]

	@test sprint(show, μ) == "μ"

    @objective gmp Max Mom(μ, x)
    cref = @constraint gmp Mom(x,μ) + Mom(x,ν) == 1
    @test sprint(show, cref) == "⟨μ, x⟩ + ⟨ν, x⟩ = 1.0"
    @test sprint(show, gmp) == """
    GMPModel:
    Maximize ⟨μ, x⟩
    s.t.
    ⟨μ, x⟩ + ⟨ν, x⟩ = 1.0
    Unknowns: μ ν
    OPTIMIZE_NOT_CALLED
    """
    @objective gmp Min Mom(μ, x)
    @constraint gmp Mom(x,μ) + Mom(x,ν) <= 1
    mc = @constraint gmp Mom(x,μ) + Mom(x,ν) >= 1
    @test sprint(show, mc) == "⟨μ, x⟩ + ⟨ν, x⟩ ≥ 1.0"
    @test sprint(show, gmp) == """
    GMPModel:
    Minimize ⟨μ, x⟩
    s.t.
    ⟨μ, x⟩ + ⟨ν, x⟩ = 1.0
    ⟨μ, x⟩ + ⟨ν, x⟩ ≤ 1.0
    ⟨μ, x⟩ + ⟨ν, x⟩ ≥ 1.0
    Unknowns: μ ν
    OPTIMIZE_NOT_CALLED
    """


end  
