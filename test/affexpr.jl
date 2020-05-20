@testset "Expressions" begin
    @testset "AbstractGMPExpr Types" begin

        @test MO.GMPEmptyExpr() isa MO.AbstractGMPExpressionLike

        @testset "ObjectExpr" begin
            @polyvar x y
            m = GMPModel()
            @variable m μ Meas([x])
            @variable m ν Meas([x])
            @variable m σ Meas([y], basis = ChebyshevBasis)
            @test_throws AssertionError MO.same_variables([μ, σ])
            @test MO.same_variables([μ, ν]) isa Nothing
            @test MO.same_vref_type([μ, σ]) isa Nothing

            oe1 = μ + ν

            @test MO.gmp_coefficients(μ) == [1]
            @test MO.gmp_coefficients(oe1) == [1, 1]
            @test MO.gmp_variables(μ) == [μ]
            @test MO.gmp_variables(oe1) == [μ, ν]

            @test !iszero(oe1)
            @test oe1[μ] == 1
            @test length(oe1) == 2
            @test collect(oe1) isa AbstractVector
            @test oe1 == oe1
            @test !(oe1 == ν + μ) #ObjectExpr do not commute

            @test variables(oe1) == [x]
            @test MO.vref_type(oe1) == MO.AbstractGMPMeasure
            @test MO.degree(oe1) == 0

            @test sprint(show, MO.ObjectExpr(Int[], MO.GMPVariableRef{MO.AbstractGMPMeasure}[])) == "0"
            o2 = (μ + ν)/2
            @test eltype([μ, o2]) == MomentOpt.ObjectExpr{typeof(1.0), MO.AbstractGMPMeasure}
            @test eltype([oe1, o2]) == MomentOpt.ObjectExpr{typeof(1.0), MO.AbstractGMPMeasure}

        end
        @testset "MomentExpr" begin
            @polyvar x y
            m = GMPModel()
            @variable m μ Meas([x])
            @variable m ν Meas([x])
            @variable m σ Meas([x]; basis = ChebyshevBasis)

            m1 = MomentExpr(1, μ)
            @test sprint(show, m1) == "⟨1, μ⟩"
            @test MO.gmp_coefficients(m1) == [1]
            @test MO.gmp_variables(m1) == [1*μ] 
            @test !(iszero(m1))
            @test m1[μ] isa Nothing
            @test m1[1*μ] == 1
            @test length(m1) == 1
            @test collect(m1) isa AbstractVector
            @test m1 == m1
            o1 = μ
            o2 = (μ+ν)/2
            m2 = MomentExpr([0.5, 0.5], [μ, ν])
            m3 = Mom(1, o2)
            @test !(m2 == m3)
            @test sprint(show, m2) == "⟨0.5, μ⟩ + ⟨0.5, ν⟩"
            @test sprint(show, m3) == "⟨1, 0.5μ + 0.5ν⟩" 
            @test sprint(show, m1 - m1) == "⟨0, 0⟩"
            @test MO.momexp_by_measure(m2) == MO.momexp_by_measure(m3)
            m4 = Mom([1, x], [o1, o2])        
            m5 = Mom([-x, 1],[o1, o2])
            @test MO.degree(m1) == 0 
            @test MO.degree(m2) == 0
            @test MO.degree(m3) == 0
            @test MO.degree(m4) == 1

            @test sprint(show, m5) == "⟨-x, μ⟩ + ⟨1, 0.5μ + 0.5ν⟩"
            @test sprint(show, m4 - m5) == "⟨x + 1, μ⟩ + ⟨x - 1, 0.5μ + 0.5ν⟩" 
            @test sprint(show, MomentExpr(Int[], MO.GMPVariableRef{MO.AbstractGMPMeasure}[])) == "⟨0, 0⟩"

            @test eltype([Mom(1, μ + ν), μ]) == AbstractJuMPScalar
            @test eltype([Mom(1, μ + ν), μ + ν]) == MO.AbstractGMPExpr
            @test eltype([Mom(1, μ + ν), Mom(one(Int8), μ + ν)]) == MomentExpr{Int, Int}
            
            @test_throws ErrorException integrate(m5)
            z = zero(m5)
            @test isempty(MO.gmp_coefficients(z))
            @test isempty(MO.gmp_variables(z))
        end

    end

    @testset "AbstractGMPExpr - Arith" begin

        @testset "ObjectExpr - Arith" begin
            @polyvar x y
            m = GMPModel()
            @variable m μ Meas([x])
            @variable m ν Meas([x])
            @variable m σ Meas([x]; basis = ChebyshevBasis)
            o1 = μ
            o2 = (μ+ν)/2
            @test sum([o1, o2])+ν == 3*o2
            @test o1-μ == Mom(Int[], MO.GMPVariableRef{MO.AbstractGMPMeasure}[])
            o3 =(μ-ν)/2
            @test o2 + o3 == MO.ObjectExpr(1, o1)
            @test 2*o2 == o2*2

            λ = lebesgue_measure_box(variable_box(x => [0, 1], y => [0, 1]))
            @test sprint(show, μ + λ) == "μ + AnalyticMeasure"
            @test sprint(show, μ + λ - λ) == "μ + AnalyticMeasure"
            @test sprint(show, 2*(μ + λ)) == "2μ + AnalyticMeasure"
            @test sprint(show, (μ + λ)/4) == "0.25μ + AnalyticMeasure"


        end

        @testset "MomentExpr - Arith" begin
            @polyvar x y
            m = GMPModel()
            @variable m μ Meas([x])
            @variable m ν Meas([x])
            @variable m σ Meas([x]; basis = ChebyshevBasis)
            o1 = μ
            o2 = (μ+ν)/2
            o3 =(μ-ν)/2

            m1 = MomentExpr(1.0, μ)        
            m2 = MomentExpr([0.5, 0.5], [μ, ν])
            m3 = Mom(1, o2)

            @test m1 + Mom(1, μ) == Mom(2, μ)
            @test length(sum([m1, m2, m3])) == 3
            @test length(m3 + m1 - 2*m2) == 2
            @test MO.momexp_by_measure(Mom(1, o2) + Mom(1, o3)) == m1

            mv = Mom.([1, 2], μ)
            @test JuMP.isequal_canonical(mv[2] - mv[1], mv[1])
            @test Mom.(1, [μ, μ]) isa Vector{MomentExpr{Int, Int}}
            ae1 = m1 + 1
            ae2 = 2*ae1 - 1
            ae3 = m1 - 1
            @test MO.constant(ae1) == MO.constant(ae2)
            @test ae2*2 - ae2 == ae2
            @test sprint(show, ae2*2 - 2) == "⟨4.0, μ⟩"
            @test ae1 - 1 == -(-ae1 + 1)
        end

    end
end
