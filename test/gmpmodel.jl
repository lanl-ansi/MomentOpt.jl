@testset "gmpmodel" begin
    @test supertype(GMPModel) == JuMP.AbstractModel
    @test GMPModel() isa GMPModel

    gmp = GMPModel()
    @test MomentOpt.model_info(gmp) isa MomentOpt.ModelInfo
    @test JuMP.num_variables(gmp) == 0
    @test isempty(JuMP.object_dictionary(gmp))
    @test MomentOpt.relax_info(gmp) isa MomentOpt.RelaxInfo
    @test MomentOpt.relax_model(gmp) isa Nothing
end
