export approximate!
"""
    approximate!(m::GMPModel)

Computes an approximation of m.
"""
function approximate!(model::GMPModel)
    model.approximation_model = Model()
    set_optimizer(approximation_model(model), approximation_info(model).solver) 
    return approximate!(model, approximation_mode(model))
end

function approximate!(model::GMPModel, mode::AbstractPrimalMode)
    # init approx objects
    # substitutions/delete JuMP variables (better when creating) (substitutions should only work with objects, that are already defined)
    # if measure 
    #   add justification constraints justify(vref, Approximationtype, mode)
    # else
    #   for continuous functions the approximation is not justified
    # end
    #
    # add constraints
    # if nonneg constraint 
    #   add justification
    # end
    # add objective
    optimize!(approximation_model(model))

    # set dual justifications 
    return nothing
end

function dual_variable(model::JuMP.Model, con::AbstractGMPConstraint, sense::MOI.OptimizationSense, deg::Int)
    if con isa MeasureConstraint
        # for measure constraints the basis of the first registered measure is used to test equality should probably be moved to the analytic measure...
        v = vref_object(first(gmp_variables(jump_function(con))))
        variable = @variable(model, variable_type = Poly(maxdegree_basis(v, deg)))
    elseif moi_set(con) isa MOI.EqualTo
        variable = @variable(model)
    elseif xor(moi_set(con) isa MOI.LessThan, sense == MOI.MIN_SENSE)
        variable = @variable(model, lower_bound=0)
    else
        variable = @variable(model, upper_bound=0)
    end
    return variable
end

function approximate!(model::GMPModel, mode::AbstractDualMode)
    # init constraint multipliers
    dvar = Dict()
    for (i, con) in gmp_constraints(model)
        dvar[i] = dual_variable(approximation_model(model), con, objective_sense(model), approximation_degree(model))
    end
    
    PT = polynomialtype(eltype(variables(gmp_variables(model)[1].v)), JuMP.AffExpr)
    dlhs = Dict{Int, PT}()

    for i in keys(gmp_variables(model))
        dlhs[i] = zero(PT)        
    end
    d_obj =0

    for (i, con) in gmp_constraints(model)
        if con isa MeasureConstraint
            for v in gmp_variables(jump_function(con))
                dlhs[index(v)] = MA.add!(dlhs[index(v)], dvar[i])
            end
            d_obj = MA.add!(d_obj, integrate(dvar[i], MOI.constant(moi_set(con))))
        elseif con isa MomentConstraint
            mcon = momexp_by_measure(jump_function(con))
            for (c, v) in mcon
                idx = index(first(gmp_variables(v)))
                dlhs[idx] = MA.add_mul!(dlhs[idx], dvar[i], c)
            end
            d_obj = MA.add_mul!(d_obj, dvar[i], MOI.constant(moi_set(con)))
        end
    end

    primal_sense = objective_sense(model)    
    if primal_sense == MOI.FEASIBILITY_SENSE
        obj = zero(MomentExpr{Int, Int})
        dual_sense = primal_sense
    elseif primal_sense == MOI.MAX_SENSE
        obj = - objective_function(model)
        dual_sense = MOI.MIN_SENSE
    elseif primal_sense == MOI.MIN_SENSE
        obj = objective_function(model)
        dual_sense = MOI.MAX_SENSE
    end
    obj_const = constant(obj)
    obj_expr = momexp_by_measure(expr(obj))
    dcon = Dict()
    just = Dict()
    for (i, v) in gmp_variables(model)
        p = dlhs[i]
        if primal_sense == MOI.MIN_SENSE
            p = -p
        end
        idx = findfirst(x -> index(first(gmp_variables(x))) == i, gmp_variables(obj_expr))
        if !(idx isa Nothing)
            p = MA.add!(p, gmp_coefficients(obj_expr)[idx])
        end
        multipliers = approximation_scheme(approx_scheme(v.v), bsa_set(v.v), variables(v.v), approximation_degree(model))
        just[i] = Dict()
        for (mons, set) in multipliers
            if size(mons, 1) == 1
                just[i][mons] = @variable(approximation_model(model), lower_bound = 0 )
                p = MA.sub_mul!(p, first(mons), just[i][mons]) 
            else
            just[i][mons] = @variable(approximation_model(model), [1:size(mons,1),1:size(mons,2)], PSD)
            p = MA.sub!(p, LinearAlgebra.tr(mons*just[i][mons])) 
        end
        end
        dcon[i] = @constraint approximation_model(model) p == 0
    end

    @objective(approximation_model(model), dual_sense,  d_obj)

    optimize!(approximation_model(model))

    for i in keys(gmp_variables(model))
        approx_vrefs(model)[i] = RefApprox(dual(dcon[i]), value(dcon[i]), just[i])
    end
    for i in keys(gmp_constraints(model))
        approx_crefs(model)[i] = RefApprox(nothing, value(dvar[i]), nothing)
    end
    return nothing
end

# inner approximations of positive polynomials

function QuadraticModuleMonomials(o::AbstractGMPObject, degree::Int)
    eqs = equalities(MOI.get(o, BSASet()))
    ineqs = [_mono_one(MOI.get(o, Variables())), inequalities(MOI.get(o, BSASet()))...]
    QM = Dict(:eq => Dict(), :ineq => Dict())
    # TODO reflect basis of Approximation
    for eq in eqs
        QM[:eq][eq] = monomials(MOI.get(o, Variables()), 0:degree-maxdegree(eq))
    end
    for ineq in ineqs 
        QM[:ineq][ineq] = monomials(MOI.get(o, Variables()), 0:Int(floor((degree-maxdegree(ineq))/2)))
    end
    return QM
end

function MomentPutinar(model::GMPModel, degree::Int)
    for (i, v_approx) in approximation_info(model).approx_vrefs
        if MOI.get(gmp_variables(model)[i].v, ApproximationType()) isa MomentOpt.DefaultApproximation
            qm = QuadraticModuleMonomials(gmp_variables(model)[i].v, degree)
            for (eq, mults) in qm[:eq]
                for mult in mults
                    cref =  @constraint approximation_model(model) riesz(model, i, eq*mult) == 0
                    push!(v_approx.cref, cref)
                end
            end
            for (ineq, mult) in qm[:ineq]
                if length(mult) == 1
                    moiset = MOI.GreaterThan(0.0)
                    mult = first(mult)
                else
                    moiset = PSDCone()
                end
                cref = @constraint approximation_model(model) riesz.(model, i, ineq*mult*transpose(mult)) in moiset
                push!(v_approx.crefs, cref)
            end
        end
    end
    return nothing
end


function generate_approximation_model(model::GMPModel)
    return measure_relaxation_model(model::GMPModel) # TODO be general
end

function measure_relaxation_model(model::GMPModel)
    degree = approximation_info(model).degree
    # initiate moments
    init_approximation(model, degree)
    # add substitutions
    # add measure condition to moments
    MomentPutinar(model, degree)

    # from now on no new variables shoud be added. 

    # add objective
    @objective approximation_model(model) objective_sense(model) riesz(model, objective_function(model)) 
    # add moment constraints
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MomentConstraintShape
            approximation_info(model).approx_crefs[i] = @constraint approximation_model(model) riesz(model, jump_function(con)) in moi_set(con)
        end
    end
    # add measure constraints
    return nothing
end

