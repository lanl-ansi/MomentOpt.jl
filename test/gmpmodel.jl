@testset "GMPModel" begin
    @testset "Model" begin
        model = GMPModel()
        @test MO.model_info(model) isa MO.ModelInfo
        @test MO.model_data(model) isa MO.ModelData
        @test MO.approximation_info(model) isa MO.ApproximationInfo
        @test MO.approximation_model(model) isa JuMP.Model
        @test MO.gmp_variables(model) isa Dict{Int, MO.AbstractGMPVariable}
        @test MO.gmp_constraints(model) isa Dict{Int, MO.AbstractGMPConstraint}
        @test MO.variable_names(model) isa Vector{AbstractString}
        @test MO.constraint_names(model) isa Vector{AbstractString}
        @test JuMP.object_dictionary(model) isa Dict{Symbol, Any}
        @test JuMP.variable_type(model) == MO.GMPVariableRef
        @test JuMP.constraint_type(model) == MO.GMPConstraintRef

        @polyvar x y
        @variable model [1:3] Meas([x,y]) 
        @test JuMP.num_variables(model) == 3

    end

    @testset "Objective" begin

        @polyvar x y
        m = GMPModel()
        mus = @variable m [1:2] Meas([x,y]) 
        JuMP.set_objective(m, MOI.MAX_SENSE, sum(Mom.(1, mus)))

        @test JuMP.objective_sense(m) == MOI.MAX_SENSE
        JuMP.set_objective_sense(m, MOI.MIN_SENSE)
        @test JuMP.objective_sense(m) == MOI.MIN_SENSE
        @test JuMP.objective_function(m) == Mom(1, mus[1]) + Mom(1, mus[2])
        @test JuMP.objective_function_type(m) == MO.MomentExpr{Int, Int} 
        @test JuMP.objective_function(m, MO.MomentExpr{Int, Int}) == Mom(1, mus[1]) + Mom(1, mus[2])
       end
    @testset "Summary" begin
        @polyvar x y
        m = GMPModel()
        mus = @variable m [1:2] Meas([x,y]) 
        JuMP.set_objective(m, MOI.MAX_SENSE, sum(Mom.(1, mus)))
        @test sprint(show, m) == "A JuMP Model\nMaximization problem with:\nVariables: 2\nObjective function type: MomentExpr{Int64,Int64}\nConstraints: 0\nApproxmation mode: DUAL_STRENGTHEN_MODE()\nMaximum degree of data: 0\nDegree for approximation 0\nSolver for approximation: "
        
        @constraint m Mom(1, mus[1]) <= 1
        @constraint m c2 Mom(1, mus[2]) >= 1
        @test sprint(print, m) == "Max ⟨1, noname⟩ + ⟨1, noname⟩\nSubject to\n ⟨1, noname⟩ ≤ 1.0\n ⟨1, noname⟩ ≥ 1.0\n"
    end

    @testset "Get Value" begin

        @polyvar x y
        m = GMPModel()
        @variable m mu Meas([x])     
        @test_throws AssertionError approximation(mu)
        @test_throws AssertionError approximation(mu, x)
        @test_throws AssertionError approximation(mu, y)
        @test_throws AssertionError JuMP.value(Mom(x, mu))
        
        @test JuMP.termination_status(m) == MOI.OPTIMIZE_NOT_CALLED

        @test MO.approximation_degree(m) == MO.model_degree(m)
        @test_logs (:warn, "Requested approximation degree -1 is too low to cover all data. The approximation degree has been set to the minimal value possible and now is 0.") set_approximation_degree(m, -1)
        @test MO.approximation_degree(m) == 0
        set_approximation_degree(m, 1)
        @test MO.approximation_degree(m) == 2

        @test MO.approximation_mode(m) isa MO.DUAL_STRENGTHEN_MODE
        set_approximation_mode(m, MO.PRIMAL_RELAXATION_MODE())
        @test MO.approximation_mode(m) isa MO.PRIMAL_RELAXATION_MODE

    end

end

