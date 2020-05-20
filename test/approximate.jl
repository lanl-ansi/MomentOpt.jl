function testproblem1(factory)
    @polyvar x y 
    f = x^2*y^2 + x^4 - y^2 + 1
    gmp = GMPModel(factory)
    mu = @variable gmp μ Meas([x, y]; support = @set(x^2-x^4>=0 && 1-x^2>=0&& x-y^2==0))
    @objective gmp Min Mom(f, μ)
    con = @constraint gmp Mom(1, μ) == 1
    set_silent(gmp)
    return gmp, mu, con
end

function testproblem2(factory)
    @polyvar x y 
    gmp = GMPModel(factory)
    mu = @variable gmp μ Meas([x, y]; support = @set(1-x^2-y^2>=0))
    @variable gmp ν Meas([x, y]; support = @set(1-x^2>=0 && 1-y^2>=0))
    @objective gmp Max Mom(1, μ)
    con = @constraint gmp μ + ν == lebesgue_measure_box( x => [-1, 1], y => [-1, 1]) 
    set_silent(gmp)
    return gmp, mu, con
end

@testset "Approximate" begin
    @testset "Optimize not called" begin
        m, mu, con = testproblem1(SCS.Optimizer)
        @test termination_status(m) == MOI.OPTIMIZE_NOT_CALLED
        @test approximation_degree(m) == 4
        @test approximation_mode(m) == DUAL_STRENGTHEN_MODE()
        set_approximation_mode(m,  PRIMAL_RELAXATION_MODE())
        @test approximation_mode(m) == PRIMAL_RELAXATION_MODE()
        @test_throws AssertionError value(mu)
        @test_throws AssertionError primal_justification(mu)
        @test_throws AssertionError dual_justification(mu)
        @test !unset_silent(m)

        m, mu, con = testproblem2(SCS.Optimizer)
        @test approximation_degree(m) == 2
        set_approximation_degree(m, 6)
        @test approximation_degree(m) == 6
    end

    
    @testset "default dual_mode" begin
        m, mu, con = testproblem1(SCS.Optimizer)
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,N} where N, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-6))
        @test dual(con) isa Vector{Float64}

        m, mu, con = testproblem2(SCS.Optimizer)
        set_approximation_degree(m, 6)
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,2}, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-6))
        @test dual(con) isa Polynomial{true, Float64}

    end
       @testset "default primal_mode" begin
        m, mu, con = testproblem1(SCS.Optimizer)
        set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,N} where N, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-6))
        @test dual(con) isa Vector{Float64}

        m, mu, con = testproblem2(SCS.Optimizer)
        set_approximation_degree(m, 6)
        set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,2}, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-6))
        @test dual(con) isa Polynomial{true, Float64}

    end
end
