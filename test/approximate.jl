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
        @test all(isapprox.(residual(con), 0.0; atol = 1e-3))
        @test dual(con) isa Vector{Float64}

        m, mu, con = testproblem2(SCS.Optimizer)
        set_approximation_degree(m, 6)
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,2}, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-3))
        @test dual(con) isa Polynomial{true, Float64}
        
        @test isapprox(objective_value(m), integrate(1, mu); atol = 1e-3)
        @test isapprox(dual_objective_value(m), integrate(1, mu); atol = 1e-3)

        @test primal_status(m) isa MOI.ResultStatusCode
        @test dual_status(m) isa MOI.ResultStatusCode
        @test raw_status(m) isa String
        @test solve_time(m) isa Number

    
    end
    @testset "default primal_mode" begin
        m, mu, con = testproblem1(SCS.Optimizer)
        set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,N} where N, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-3))
        @test dual(con) isa Vector{Float64}

        m, mu, con = testproblem2(SCS.Optimizer)
        set_approximation_degree(m, 6)
        set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
        optimize!(m)

        @test value(mu) isa MultivariateMoments.Measure
        @test primal_justification(mu) isa Array{Array{Float64,2}, 1}
        @test dual_justification(mu) isa Polynomial{true, Float64}
        @test all(isapprox.(residual(con), 0.0; atol = 1e-3))
        @test dual(con) isa Polynomial{true, Float64}

    end

    @testset "Schmuedgen" begin
        @testset "Dual Schmuedgen" begin
            @polyvar x y 
            f = x^2*y^2 + x^4 - y^2 + 1
            gmp = GMPModel(SCS.Optimizer)
            K = @set(x^2-x^4>=0 && 1-x>=0 && x>=0 && x-y^2==0)
            mu = @variable gmp μ Meas([x, y]; support = K, scheme = SchmuedgenScheme())
            @objective gmp Min Mom(f, μ)
            con = @constraint gmp Mom(1, μ) == 1
            set_silent(gmp)

            optimize!(gmp)
            
            @test length(primal_justification(mu)) == 6
            @test size.(primal_justification(mu)) == [(6, 6), (1,), (3, 3), (3, 3),
                                                      (3, 3), (6,)] 
        end
     @testset "Primal Schmuedgen" begin
            @polyvar x y 
            f = x^2*y^2 + x^4 - y^2 + 1
            gmp = GMPModel(SCS.Optimizer)
            K = @set(x^2-x^4>=0 && 1-x>=0 && x>=0 && x-y^2==0)
            mu = @variable gmp μ Meas([x, y]; support = K, scheme = SchmuedgenScheme())
            @objective gmp Min Mom(f, μ)
            con = @constraint gmp Mom(1, μ) == 1
            set_silent(gmp)
            set_approximation_mode(gmp, PRIMAL_RELAXATION_MODE())
            optimize!(gmp)
            
            @test length(primal_justification(mu)) == 6
            @test size.(primal_justification(mu)) == [(6, 6), (1,), (3, 3), (3, 3),
                                                      (3, 3), (6,)] 
        end

    end

end


@testset "Post Proc" begin 
    @polyvar x y 
    f = x^2*y^2 + x^4 - y^2 + 1
    gmp = GMPModel()
    mu = @variable gmp μ Meas([x, y]; support = @set(x^2-x^4>=0 && 1-x^2>=0&& x-y^2==0))
    MO.approximation_info(gmp).approx_vrefs[1] = MO.VarApprox(MultivariateMoments.Measure([0.04301405019233425, 3.725430564600456e-17, 0.0944511405218325, -3.667541691651414e-17, 0.20739824193451498, 0.09445115553035512, -6.278937572688163e-17, 0.207398112278084, 2.2189511821478325e-17, 0.20739824605238844, 6.290839177153608e-18, 0.4554099572402294, 0.45540988993201603, 4.3223432627485784e-17, 1.0000002747597936], Monomial{true}[x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1]), nothing, nothing)
    @test length(atomic(mu)) == 2
    @test integrate(f, mu) == 0.682055508233731

    g =  graph(mu, [x]; regpar = 1e-6)
    @test g isa Function
    @test g([1]) isa Number
end
