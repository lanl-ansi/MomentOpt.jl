@testset "info" begin
@testset "model_info" begin
    @test supertype(MomentOpt.AbstractMeasureRef) == JuMP.AbstractVariableRef
    @test MomentOpt.NotRelaxed() isa MomentOpt.AbstractRelaxationStatus
    @test MomentOpt.PrimalRelaxed() isa MomentOpt.AbstractRelaxationStatus
    @test MomentOpt.DualRelaxed() isa MomentOpt.AbstractRelaxationStatus
    @test MomentOpt.ModelInfo() isa MomentOpt.ModelInfo
    info = MomentOpt.ModelInfo()
    @test MOI.add_variable(info) == 1
    @test MOI.add_constraint(info) == 1
end

@testset "relax_info" begin
    @test MomentOpt.PrimalFormulation() isa MomentOpt.AbstractRelaxationMode
    @test MomentOpt.DualFormulation() isa MomentOpt.AbstractRelaxationMode
    @test MomentOpt.RelaxInfo() isa MomentOpt.RelaxInfo
end
end
