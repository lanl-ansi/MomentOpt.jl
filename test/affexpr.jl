@testset "Expressions" begin
    @testset "AbstractGMPExpr Types" begin

        @testset "ObjectExpr" begin
            @polyvar x y
            m = GMPModel()
            @variable m μ Meas([x])
            @variable m ν Meas([x])
            @variable m σ Meas([y], basis = ChebyshevBasis)
            @test_throws AssertionError MO.same_variables([μ, σ])
            @test MO.same_variables([μ, ν]) isa Nothing

            oe1 = μ + ν

            @test MO.coefficients(μ) == [1]
            @test MO.coefficients(oe1) == [1, 1]
            @test MO.measures(μ) == [μ]
            @test MO.measures(oe1) == [μ, ν]

            @test !iszero(oe1)
            @test oe1[μ] == 1
            @test length(oe1) == 2
            @test collect(oe1) isa AbstractVector
            @test oe1 == oe1
            @test !(oe1 == ν + μ) #ObjectExpr do not commute

            foo = μ + μ
            @test foo[μ] == 2

            @test variables(oe1) == [x]
            @test MO.degree(oe1) == 0

            @test sprint(show, MO.MeasExpr(Int[], MO.GMPVariableRef[])) == "0"
            o2 = (μ + ν)/2
            @test_broken eltype([μ, o2]) == MomentOpt.ObjectExpr{typeof(1.0), MO.AbstractGMPMeasure}
            @test eltype([oe1, o2]) == MomentOpt.MeasExpr{typeof(1.0)}

        end
        @testset "MomentExpr" begin
            @polyvar x y
            m = GMPModel()
            @variable m μ Meas([x])
            @variable m ν Meas([x])
            @variable m σ Meas([x]; basis = ChebyshevBasis)

            m1 = MomentExpr(1, μ)
            @test sprint(show, m1) == "⟨1, μ⟩"
            @test MO.coefficients(m1) == [1]
            @test MO.measures(m1) == [μ] 
            @test !(iszero(m1))
            @test m1[μ] == 1
            @test length(m1) == 1
            @test collect(m1) isa AbstractVector
            @test m1 == m1
            o1 = μ
            o2 = (μ+ν)/2
            m2 = MomentExpr([0.5, 0.5], [μ, ν])
            m3 = Mom(1, o2)
            @test (m2 == m3)
            @test sprint(show, m2) == "⟨0.5, μ + ν⟩"
            @test sprint(show, m3) == "⟨0.5, μ + ν⟩"
            @test sprint(show, m1 - m1) == "⟨0, 0⟩"
            m4 = Mom([1, x], [o1, o2])        
            m5 = Mom([-x, 1],[o1, o2])
            @test MO.degree(m1) == 0 
            @test MO.degree(m2) == 0
            @test MO.degree(m3) == 0
            @test MO.degree(m4) == 1

            @test sprint(show, m5) == "⟨-x + 0.5, μ⟩ + ⟨0.5, ν⟩"
            @test sprint(show, m4 - m5) == "⟨1.5x + 0.5, μ⟩ + ⟨0.5x - 0.5, ν⟩"
            @test sprint(show, MomentExpr(Int[], MO.GMPVariableRef[])) == "⟨0, 0⟩"

            @test eltype([Mom(1, μ + ν), μ]) == AbstractJuMPScalar
            @test eltype([Mom(1, μ + ν), μ + ν]) == MO.AbstractGMPScalar
            @test eltype([Mom(1, μ + ν), Mom(one(Int8), μ + ν)]) == MomentExpr{Int}
            
            z = zero(m5)
            @test isempty(MO.coefficients(z))
            @test isempty(MO.measures(z))
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
            @test o1-μ == MO.MeasExpr(Int[], MO.GMPVariableRef[])
            o3 =(μ-ν)/2
            @test o2 + o3 == MO.MeasExpr(1, o1)
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
            @test length(sum([m1, m2, m3])) == 2 
            @test length(m3 + m1 - 2*m2) == 2

            mv = Mom.([1, 2], μ)
            @test_broken JuMP.isequal_canonical(mv[2] - mv[1], mv[1])
            @test Mom.(1, [μ, μ]) isa Vector{MomentExpr{Int}}
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
