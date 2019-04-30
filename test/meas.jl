@testset "Measure Tests" begin
	# test  Measure, variables, support, certificate
	@polyvar x y
	K = @set(x^2+y^2<=1 && x*y >=0)
	μ = Measure("μ", [x,y],support = K)
	@test variables(μ) == [x,y]
	@test support(μ) == K
	@test certificate(μ) == SOSCone()
	ν = Measure("ν",[x],certificate=DSOSCone())
	@test support(ν) ==  FullSpace()
	@test certificate(ν) == DSOSCone()
end
