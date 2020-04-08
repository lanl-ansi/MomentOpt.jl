@testset "Objects" begin
    @testset "AbstractGMPObject" begin

        @polyvar x y

        mone = MO._mono_one([x, y])
        @test isone(mone)
        @test mone isa Monomial{true}
        @test variables(mone) == [x, y]

        mone = MO._mono_one(nothing) 
        @test mone isa Number
        @test isone(mone)

    end

    @testset "AbstractGMPMeasure" begin

        @polyvar x y
        p1 = x^2 + 4*x + 1
        p2 = x*y + 2*(x + y) +1
        S1 = @set x^2 <= 1
        S2 = @set x^2 <= 1 && y == x

        z = ZeroMeasure()
        @test sprint(show, z) == "0" 
        @test monomials(MO.covering_basis(z, p1)) == [x^2, x, 1] 
        @test integrate(1, z) == 0
        @test integrate(p1, z) == 0

        #AbstractGMPObjects should check whether support and variables are coherent

        σ = SymbolicMeasure([x], S1, MonomialBasis)
        @test MO.covering_basis(z, p1) == MO.covering_basis(σ, p1)
        @test !(z==σ)
        @test_throws AssertionError integrate(p1, σ)

        λ = AnalyticMeasure([x], MonomialBasis, p-> sum(c/(1+first(exponents(m))) for (c,m) in zip(coefficients(p), monomials(p))))
        @test λ == λ
        @test integrate(p1, λ) == 10/3
        @test_throws AssertionError integrate(p2, λ)

        μ = Meas([x, y]; support = S2, basis = ChebyshevBasis)
        @test_throws AssertionError integrate(p2, μ)
        @test monomials(maxdegree_basis(μ, 2)) == [2.0*x^2 - 1.0, x*y, 2.0*y^2 - 1.0, x, y, 1.0]  

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

    λ = lebesgue_measure_box(box, basis = MonomialBasis)
    @test λ isa AnalyticMeasure
end
end
