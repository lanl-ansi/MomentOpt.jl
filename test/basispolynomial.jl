@testset "basispolynomial" begin
    @polyvar x y
    f = x^3 + x^2*y - x^2 - x*y + x - 1

    fb = basis_polynomial(MonomialBasis, f)
    @test coefficients(fb) == [1, 1, -1, -1, 1, -1]
    @test monomials(fb) == monomials(f)
    @test polynomial(fb) == f

    fb = basis_polynomial(LaguerreBasis, f)
    @test coefficients(fb) == [-6, -2, 18, 3, -18, -1, 5]
    @test length(monomials(fb)) == 7
    @test polynomial(fb) == f

    fb = basis_polynomial(MonomialBasis, y + x^2)
    @test sprint(show, fb) == "1.0 * ( x^2 ) + 1.0 * ( y )"
end

@testset "riesz functional" begin
    @polyvar x y
    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x], 1), ones(2))
    @test sprint(show, ms) == "Dict(x=>1.0,1.0=>1.0)"
    @test_throws AssertionError moment_sequence(maxdegree_basis(ChebyshevBasis, [x], 2), ones(4))

    f = x^3 + x^2*y - x^2 - x*y + x - 1
    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x, y], 3), ones(10))
    @test riesz(ms, f) == sum(coefficients(f))
    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x, y], 2), ones(6))
    @test_throws AssertionError riesz(ms,f)
    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x], 2), [1, 0, 1])
    @test_throws AssertionError riesz(ms,f)
    mv = monomials([x], 0:1)
    mm = mv*transpose(mv)
    @test riesz.(ms, mm) == [1.0 0.0; 0.0 1.0]

    @polyvar z[1:3]
    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x], 2), [z[3], z[2], z[1]])
    mm = mv*transpose(mv)

    @test riesz.(ms, mm) == [0.5*(z[1]+z[3]) z[2]; z[2] z[1]]

    ms = moment_sequence(maxdegree_basis(ChebyshevBasis, [x, y], 3), monomials(ChebyshevBasis, [x, y], 0:3))
    @test riesz(ms, f) == f
    ms = moment_sequence(maxdegree_basis(MonomialBasis, [x], 4), [1/5, 1/4, 1/3, 1/2, 1])
    mv = monomials(maxdegree_basis(LaguerreBasis, [x], 2))
    mm = mv*transpose(mv)
    @test riesz.(ms, mm) ≈ [13/60 25/120 1/6; 25/120 1/3 1/2; 1/6 1/2 1]

end

