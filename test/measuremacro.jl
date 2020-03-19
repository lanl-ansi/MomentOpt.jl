@testset measure_macro begin
    m = GMPModel()
    @polyvar x y

    @test_throws ArgumentError @measure m
    
    unnamed = @measure m [x, y]
    @test unnamed isa MeasureRef
    @test JuMP.name(unnamed) == ""
    @test poly_variables(unnamed) == [x, y]
    @test support(unnamed) == FullSpace()
    @test moment_basis(unnamed) == MonomialBasis
    @test relax_type(unnamed) isa MomentOpt.Putinar #TODO: export Putinar 

    @measure m µ [x]
    @test owner_model(µ) == m
    @test variable_by_name(m, "μ") == µ
    @test poly_variables(μ) == [x]
    
    set_poly_variables(μ, [x,y])
    @test poly_variables(μ) == [x, y]
    set_support(μ, @set( 1-x^2-y^2 >=0))
    @test support(μ) isa BasicSemialgebraicSet
    set_moment_basis(μ, ScaledMonomialBasis)
    @test moment_basis(μ) == ScaledMonomialBasis
    set_relax_type(μ, MomentOpt.Handelman()) #TODO: export Handelman
    @test relax_type(μ) isa MomentOpt.Handelman #TODO: export Handelman

    ν = @measure m [1:2] [x, y] 
    @test length(ν) == 2
    @test support(ν[1]) == support(ν[2])
    @test set_relax_type.(ν, MomentOpt.Handelman) isa Vector{MeasureRef}  #TODO: export Handelman
    
    @test @measure(m, ϕ[1:2], [x]) isa Vector{MeasureRef}
    
    # keywords
    K = @set( 1-x^2-y^2 >=0 )
    foo = @measure m [x,y] support = K
    @test support(foo) == K

    foo = @measure m [x,y] moment_basis = ChebyshevBasis
    @test moment_basis(foo) == ChebyshevBasisFirstKind

    foo = @measure m [x,y] relax_type = MomentOpt.Handelman #TODO: export Handelman
    @test relax_type(foo) == MomentOpt.Handelman #TODO: export Handelman

    foo = @measure m [x,y] support = K moment_basis = ChebyshevBasis relax_type = MomentOpt.Handelman #TODO: export Handelman
    @test support(foo) == K
    @test moment_basis(foo) == ChebyshevBasisFirstKind
    @test relax_type(foo) == MomentOpt.Handelman #TODO: export Handelman

    @measure m [1:2] [x] base_name = "bar"
    @test bar isa Vector{MeasureRef}

end
