
function JuMP.optimize!(model::GMPModel)
    # TODO instead of creating a new model we should rather call empty!
    # We can than remove the solver field of approximation_info and make set_optimizer/ optimizer_with_attributes, set_optimizer_attribute(s), get_optimizer_attribute, pass directly to the jumo model. 
    model.approximation_model = Model(approximation_info(model).solver)
    return approximate!(model, approximation_mode(model))
end

function approximate!(model::GMPModel, mode::AbstractPrimalMode)
    #TODO clean up
    approximate_putinar!(model, mode)
    return nothing
end

function dual_variable(model::JuMP.Model, con::AbstractGMPConstraint, sense::MOI.OptimizationSense, deg::Int)
    if con isa MeasureConstraint
        # for measure constraints the basis of the first registered measure is used to test equality. it should probably be moved to the analytic measure in the moiset
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
    # init dual variables
    dvar = Dict()
    for (i, con) in gmp_constraints(model)
        dvar[i] = dual_variable(approximation_model(model), con, objective_sense(model), approximation_degree(model))
    end
    
    # init dual constraints
    PT = polynomialtype(eltype(variables(gmp_variables(model)[1].v)), JuMP.AffExpr)
    dlhs = Dict{Int, PT}()
    for i in keys(gmp_variables(model))
        dlhs[i] = zero(PT)        
    end
    d_obj =0

    # summing over measures for each constraint. 
    for (i, con) in gmp_constraints(model)
        if con isa MeasureConstraint
            for v in gmp_variables(jump_function(con))
                dlhs[index(v)] = MA.add!(dlhs[index(v)], dvar[i])
            end
            d_obj = MA.add!(d_obj, integrate(dvar[i], moi_set(con)))
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
        scheme_parts = approximation_scheme(model, v)
        just[i] = [dual_scheme_variable(approximation_model(model), sp)*sp.polynomial for sp in scheme_parts]
        p = MA.sub!(p, sum(just[i] )) 
        dcon[i] = @constraint approximation_model(model) p == 0
    end

    # set objective
    @objective(approximation_model(model), dual_sense, d_obj)

    # call optimize on the approximation model
    optimize!(approximation_model(model))

    # store values and duals
    for i in keys(gmp_variables(model))
        moms = dual(dcon[i])
        approx_vrefs(model)[i] = RefApprox(moms, value(dcon[i]), value.(integrate.(just[i], moms)))
    end
    for (i, con) in gmp_constraints(model)
        if dvar[i] isa MP.AbstractPolynomialLike
            pvalue =  0 #integrate.(dvar[i], jump_function(con))
        elseif has_lower_bound(dvar[i])
            pvalue = dual(LowerBoundRef(dvar[i]))
        elseif has_upper_bound(dvar[i])
            pvalue = dual(UpperBoundRef(dvar[i]))
        else
            pvalue = 0
        end
        approx_crefs(model)[i] = RefApprox(pvalue, value(dvar[i]), nothing)
    end
    return nothing
end

function moments_variable(model::JuMP.Model, v::GMPVariable{AbstractGMPMeasure}, deg::Int)
    X = maxdegree_basis(approx_basis(v.v), variables(v.v), deg)
    c = @variable model [1:length(X)] 
    return measure(c, monomials(X))
end

function approximate_putinar!(model::GMPModel, ::AbstractPrimalMode)
    degree = approximation_info(model).degree
    # initiate moments 
    pvar = Dict(i => moments_variable(approximation_model(model), v, degree) for (i,v) in gmp_variables(model))
    # add substitutions

    # add measure condition on moments
    just =  Dict{Int, Vector{ConstraintRef}}()
    for (i, v) in gmp_variables(model)
        just[i] = ConstraintRef[]
        scheme_parts = approximation_scheme(model, v)
        for sp in scheme_parts
            cref = primal_scheme_constraint(approximation_model(model), sp, pvar[i])
            if cref isa AbstractVector
                append!(just[i], cref)
            else
                push!(just[i], cref)
            end
        end
    end

    # add constraints
    pcon = Dict{Int, Vector{ConstraintRef}}()
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MomentConstraintShape
            cref = @constraint approximation_model(model) integrate(jump_function(con), pvar) in moi_set(con)            
            pcon[i] = [cref]
        elseif shape(con) isa MeasureConstraintShape
        # add measure constraints
        refmeas = moi_set(con)
        mons = monomials(maxdegree_basis(approx_basis(refmeas), variables(refmeas), approximation_degree(model)))
        pcon[i] = @constraint approximation_model(model) sum(c.*(integrate.(mons, pvar[index(v)])) for (c,v) in jump_function(con)) .== integrate.(mons, refmeas) 
        end
    end
       # add objective
    @objective approximation_model(model) objective_sense(model) integrate(objective_function(model), pvar)

    optimize!(approximation_model(model))

   for i in keys(gmp_variables(model))
       approx_vrefs(model)[i] = RefApprox(MM.Measure(value.(pvar[i].a), pvar[i].x), nothing, value.(just[i]))
    end
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MeasureConstraintShape
        refmeas = moi_set(con)
        mons = monomials(maxdegree_basis(approx_basis(refmeas), variables(refmeas), approximation_degree(model)))

        #todo check whether minus has something to do with objective sense
        approx_crefs(model)[i] = RefApprox(value.(pcon[i]), -dot(dual.(pcon[i]), mons), nothing)
        else
            approx_crefs(model)[i] = RefApprox(value.(pcon[i]), dual.(pcon[i]), nothing)
        end
    end

    return nothing
end

