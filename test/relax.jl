@testset "Relax Test" begin
    @polyvar x y
    f =x^4*y^2 + x^2*y^4 -3x^2*y^2 + 1 
    K = @set(1-x>=0 && x+1>=0 && 1-y>=0 && y+1>=0)
    gmp = GMPModel()
    @measure gmp μ [x,y] support=K
    @objective gmp Min Mom(f,μ)
    @constraint gmp c Mom(1,μ) == 1
    relax!(gmp, 2, with_optimizer(CSDP.Optimizer, printlevel=0))
    @test value(gmp, Mom(μ,1)) isa Float64
    @test dual_value(gmp, c) isa Float64
    @test atomic(gmp, μ, tol = 1e-03, print_level = 0) == nothing
    @test atomic(gmp, μ, tol = 1e-03, print_level = 1) == nothing
    ch =  christoffel(gmp, μ)
    @test ch isa Polynomial{true,Float64}

    ERR = ErrorException("")
    try min_val(([x,y]=>[0,0]), ch)
    catch err
        @test err==ERR
    end

    @test min_val(([x]=>[0]), ch) isa Float64

    @test objective_value(gmp) isa Float64
    @test value(gmp,Mom(μ,x^2+y+1)) isa Float64
    relax!(gmp, 3, with_optimizer(CSDP.Optimizer, printlevel=0))
    @test atomic(gmp, μ, tol = 1e-03, print_level = 0) isa Dict{Int, Vector{Float64}}
    @test atomic(gmp, μ, tol = 1e-03, print_level = 1) isa Dict{Int, Vector{Float64}}

    gmp = GMPModel()
    @measure gmp μ [x,y] support=K

    ERR = ErrorException("Please define an objective")
    try relax!(gmp, 2, with_optimizer(CSDP.Optimizer, printlevel=0))
    catch err
        @test err == ERR
    end
    @objective gmp Max Mom(f,μ)
    
    ERR = ErrorException( "Define at least one moment constraint!")
    try relax!(gmp, 2, with_optimizer(CSDP.Optimizer, printlevel=0))
    catch err
        @test err == ERR
    end
    
    @constraint gmp Mom(1,μ) <= 1
    @constraint gmp Mom(1,μ) >= 1
    relax!(gmp, 2, with_optimizer(CSDP.Optimizer, printlevel=0))





end

