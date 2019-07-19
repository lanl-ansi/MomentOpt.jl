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
    @test atomic(gmp, μ, 1e-03) == nothing
    @test objective_value(gmp) isa Float64
    @test value(gmp,Mom(μ,x^2+y+1)) isa Float64
    relax!(gmp, 3, with_optimizer(CSDP.Optimizer, printlevel=0))
    @test atomic(gmp, μ, 1e-03) isa Dict{Int, Vector{Float64}}
end
