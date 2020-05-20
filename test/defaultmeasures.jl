function test_moments(μ::AnalyticMeasure, maxdegree, moments)
    mons = monomials(maxdegree_basis(μ, maxdegree))
    for (m, mm) in zip(mons, moments)
        @test isapprox(integrate(m, μ), mm; atol = 1e-6)
    end
end

@testset "DefMeasures" begin

    @testset "ZeroMeasure" begin
        @polyvar x y
        ζ = ZeroMeasure([x, y])
        test_moments(ζ, 2, [0, 0, 0, 0, 0, 0])
    end

    @testset "DiracMeasure" begin
        @polyvar x
        δ1 = DiracMeasure([x], [1])
        test_moments(δ1, 4, [1, 1, 1, 1, 1])

        δ0 = DiracMeasure([x], [0])
        test_moments(δ0, 4, [0, 0, 0, 0, 1])

        δ = (δ0 + δ1)/2
        test_moments(δ, 4, [1/2, 1/2, 1/2, 1/2, 1])

        @polyvar y
        test_moments(DiracMeasure([x, y], [1, 0]), 2, [1, 0, 0, 1, 0, 1])
    end

    @testset "Lebesgue Box" begin
        @polyvar x y
        λ1 = lebesgue_measure_box(x => [0, 1], y => [0,1])
        test_moments(λ1, 2, [1/3, 1/4, 1/3, 1/2, 1/2, 1])

        λ2 = lebesgue_measure_box(x => [-1, 1], y => [0,1], normalize = true)
        test_moments(λ2, 2, [1/3, 0, 1/3, 0, 1/2, 1])
        
        test_moments(2*λ1-λ2, 2, [1/3, 1/2, 1/3, 1, 1/2, 1])
    end


    @testset "MomentMeasure" begin
        @polyvar x y
        XX = range(-1, stop = 1, length = 100)
        mons = monomials([x, y], 0:2);
        z = [1/100*(sum(XX[i]^e[1]*(2*XX[i]-1)^e[2] for i = 1:100)) for e in exponents.(mons)];
        μ = moment_measure(MultivariateMoments.Measure(z, mons))
        test_moments(μ, 2, z)
    end

    @testset "(Partial) Integration" begin
        @polyvar x y
        λ = lebesgue_measure_box(x => [0, 1])
        partial_integrate(x^2+x*y-y+1, λ) == 1/3 - 1/2*y + 1
    end

    @testset "Marginal/Product" begin
        @polyvar x y        
        λ = lebesgue_measure_box(x => [-1, 1], y => [0,1], normalize = true)

        λx = marginal(λ, [x])
        λy = marginal(λ, [y])

        mons = monomials(maxdegree_basis(λ, 2))

        μ = prod_meas(λx, λy)
        for m in mons
            @test integrate(m, λ) == integrate(m, μ)
        end 
    end
end

