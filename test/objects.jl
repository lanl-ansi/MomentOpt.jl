@testset "Objects" begin
    @testset "AbstractGMPObject" begin

        @polyvar x y

        mone = MO._mono_one([x, y])
        @test isone(mone)
        @test mone isa Monomial{true}
        @test variables(mone) == [x, y]

    end

    @testset "AbstractGMPMeasure" begin

        @polyvar x y
        p1 = x^2 + 4*x + 1
        p2 = x*y + 2*(x + y) +1
        S1 = @set x^2 <= 1
        S2 = @set x^2 <= 1 && y == x

        z = ZeroMeasure([x])
        @test monomials(MO.covering_basis(z, p1)) == [x^2, x, 1] 
        @test integrate(1, z) == 0
        @test integrate(p1, z) == 0

        #AbstractGMPObjects should check whether support and variables are coherent
        λ = AnalyticMeasure([x], MonomialBasis, p-> sum(c/(1+first(exponents(m))) for (c,m) in zip(coefficients(p), monomials(p))))
        @test λ == λ
        @test integrate(p1, λ) == 10/3
        @test_throws AssertionError integrate(p2, λ)
        @test integrate(1, λ) == 1
        μ = Meas([x, y]; support = S2, basis = ChebyshevBasis)
        @test_throws MethodError integrate(p2, μ)
        @test monomials(maxdegree_basis(μ, 2)) == [2.0*x^2 - 1.0, x*y, 2.0*y^2 - 1.0, x, y, 1.0]
        @test monomials(MO.covering_basis(μ, 2)) == Polynomial{true, typeof(1.0)}[1.0]
    end

    @testset "AbstractGMPContinuous" begin
        # Base.:(==)(m1::AbstractGMPObject, m2::AbstractGMPObject)
        # MO.covering_basis(t::AbstractGMPObject, p::MP.AbstractPolynomialLike)

        #  approximate(f::AbstractGMPContinuous, max_degree::Int)
        # SymbolicContinuous(variables::Vector{V}, domain::S, monom_basis) where {S <: AbstractBasicSemialgebraicSet, V <: MP.AbstractVariable} 
        # AnalyticContinuous(variables::Vector{V}, monom_basis, monom_function::Function) where {V <: MP.AbstractVariable} 
        # ConstantContinuous(a::Number)
        # ZeroContinuous()
        # OneContinuous()
        # approximate(f::ConstantContinuous, ::Int)
        # Cont(vars; domain = FullSpace(), basis = MonomialBasis, approx = DefaultApproximation())
        # MB.maxdegree_basis(t::AbstractGMPObject, d::Int)
        # MO.eval_vector(basis::MB.AbstractPolynomialBasis, m::AbstractGMPObject)

    end

    @testset "DefaultMeasures" begin

        @polyvar x y

        @test MO.lebesgue_line(-1, 1, 0) == 2
        @test MO.lebesgue_line(0, 1, 1) == 1/2

        box = variable_box( x => [0, 1], y => [-1/2, 1/2])
        @test MO.lebesgue_box(box, x*y) == 0
        @test MO.lebesgue_box(box, x^2) == 1/3

        box = variable_box( x => [-1, 1], y => [-1/2, 1/2])
        @test MO.lebesgue_box(box, x^2*y^2; normalize = true) == 1/36
        @test MO.lebesgue_box(box, x^2*y^2) == 1/18

        λ = lebesgue_measure_box(box, basis = MonomialBasis)
        @test λ isa AnalyticMeasure

        @test 2*λ isa AnalyticMeasure

        @test 2*integrate(x, λ) == integrate(x, λ*2) == integrate(x, λ+λ)

        ζ = ZeroMeasure([x])
        @test_throws AssertionError λ + ζ
        σ = sum([DiracMeasure([x], [i]) for i in 1:4])
        @test integrate.([x, x^2], σ) == [10, 30] 

        @test integrate(x, DiracMeasure([x], [1])/2) == 0.5
    end
end
