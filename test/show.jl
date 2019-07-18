@testset "Show Test" begin
    gmp = GMPModel()
    @test sprint(show, gmp) == """
    GMPModel:
    Feasability problem:
    s.t.
    MomCon[]
    Unknowns:
    OPTIMIZE_NOT_CALLED
    """
    @polyvar x y
    @measure gmp μ [x,y]
    @measure gmp ν [x]
    @constraint gmp Mom(x,μ)+Mom(x,ν) == 1
    @test sprint(show, gmp) == """
    GMPModel:
    Feasability problem:
    s.t.
    MomCon[⟨μ, x⟩ + ⟨ν, x⟩ = 1.0]
    Unknowns: μ ν
    OPTIMIZE_NOT_CALLED
    """
end  
