@testset "Show Test" begin
    gmp = GMPModel()
    @polyvar x y
    @test sprint(show, gmp) == """
GMPModel:
Feasability problem:
s.t.

Unknowns:
OPTIMIZE_NOT_CALLED
"""
end
