@testset "Variables" begin
    @testset "GMPVariable" begin

        @polyvar x y
        mu = Meas([x, y])
        v = GMPVariable(mu)
        @test MO.object_type(v) == VariableMeasure{PolyVar{true}, FullSpace, MonomialBasis, PutinarScheme}

        _error = _ -> ""
        info = VariableInfo(true, 1, true, 1, true, 1, true, 1, false, false)
        @test JuMP.build_variable(_error, info, mu) isa MO.AbstractGMPVariable

    end

    @testset "GMPVariableRef" begin
        m = GMPModel()
        vref1 = MO.GMPVariableRef(m, 1)
        vref2 = MO.GMPVariableRef(m, 2)

        @test !iszero(vref1)
        @test vref1 == vref1
        @test !(vref1 == vref2)
        @test JuMP.isequal_canonical(vref1, vref2) == (vref1 == vref2)
        @test index.([vref1, vref2]) == [1, 2]

    end

    @testset "Add/DeleteVariable" begin
        @polyvar x y
        m = GMPModel()
        @variable m μ Meas([x,y])
        nu = @variable m ν Meas([x])

        @test JuMP.name(μ) == "μ"
        @test JuMP.name(nu) == "ν"
        
        @test JuMP.variable_by_name(m, "ν") == nu
        JuMP.set_name(μ, "ν")
        @test JuMP.variable_by_name(m, "μ") isa Nothing
        @test_throws ErrorException JuMP.variable_by_name(m, "ν")

        @test JuMP.is_valid(m, μ)
        JuMP.delete(m, μ)
        @test JuMP.variable_by_name(m, "ν") == nu
    end

    @testset "Variables" begin
        @polyvar x y 
        m = GMPModel()
        mu = Meas([x,y])
        @variable m μ mu 
        @test MO.vref_object(μ) === mu
        @test variables(μ) == [x, y]
        @test support(μ) == FullSpace()
        @test domain(μ) == FullSpace()
        @test approx_basis(μ) == MonomialBasis
    end
end
