@testset "moment sequence" begin
    @polyvar x y
    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x], 1), ones(2))
    @test sprint(show, ms) == "Dict(x => 1.0,1.0 => 1.0)"
    @test_throws AssertionError moment_sequence(maxdegree_basis(ScaledMonomialBasis, [x], 2), ones(4))

    ls = lebesgue_box_sequence(variable_box(x => [-1,1], y => [-1,1]), 4; basis_type = LegendreBasis, normalize = true)
    ms = [moments(ls)[mon] for mon in monomials(ls)]
    @test last(ms)==1.0   
    @test all(ms[1:end-1] .== 0)
    
    mv = collect(maxdegree_basis(LegendreBasis, [x], 2))
    mm = mv*transpose(mv)
    @test lebesgue_box.(variable_box(x => [-1,1]), mm; normalize = true)  ≈ [1/5 0 0; 0 1/3 0; 0 0 1]
    ls = lebesgue_box_sequence(variable_box(x => [0,1]), 4)
    @test [moments(ls)[mon] for mon in monomials(ls)] == collect(5:-1:1).^(-1)

end




