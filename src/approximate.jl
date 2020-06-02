function JuMP.optimize!(model::GMPModel) 
    empty!(approximation_model(model))
    return approximate!(model, approximation_mode(model))
end

function dual_variable(model::JuMP.Model, con::AbstractGMPConstraint, sense::MOI.OptimizationSense, deg::Int)
    if con isa MeasureConstraint
        variable = @variable(model, variable_type = Poly(maxdegree_basis(moi_set(con), deg)))
    elseif moi_set(con) isa MOI.EqualTo
        variable = @variable(model)
    elseif xor(moi_set(con) isa MOI.LessThan, sense == MOI.MIN_SENSE)
        variable = @variable(model, lower_bound=0)
    else
        variable = @variable(model, upper_bound=0)
    end
    return variable
end

function approximate!(model::GMPModel, ::AbstractDualMode)
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
            for v in measures(jump_function(con))
                dlhs[index(v)] = MA.add!(dlhs[index(v)], dvar[i])
            end
            d_obj = MA.add!(d_obj, integrate(dvar[i], moi_set(con)))
        elseif con isa MomentConstraint
            mcon = jump_function(con)
            for (c, v) in mcon
                idx = index(v)
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
    obj_expr = expr(obj)
    dcon = Dict()
    just = Dict()
    scheme_parts = Dict()
    for (i, v) in gmp_variables(model)
        p = dlhs[i]
        if primal_sense == MOI.MIN_SENSE
            p = -p
        end
        ids = findall(x -> index(x) == i, measures(obj_expr))
        for id in ids
            p = MA.add!(p, coefficients(obj_expr)[id])
        end
        scheme_parts[i] = approximation_scheme(model, v)
        just[i] = [dual_scheme_variable(approximation_model(model), sp)*sp.polynomial for sp in scheme_parts[i]]

        p = MA.sub!(p, sum(just[i] )) 
        dcon[i] = @constraint approximation_model(model) p == 0
    end

    # set objective
    @objective(approximation_model(model), dual_sense, d_obj)

    # call optimize on the approximation model
    optimize!(approximation_model(model))

    # store values and duals
    for i in keys(gmp_variables(model))

        approx_vrefs(model)[i] = VarApprox(
                                dual(dcon[i]), 
                                convert.(Array{Float64}, 
                                         primal_scheme_part.(scheme_parts[i], dual(dcon[i]))), 
                                sum(value.(just[i])))
    end
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MomentConstraintShape
            approx_crefs(model)[i] = ConApprox(
                                     [integrate(Dict(i => dual(dcon[i]) for i in keys(gmp_variables(model))), 
                                                jump_function(con)) - MOI.constant(moi_set(con))], 
                                     [value.(dvar[i])])
        else
            dpol = value.(dvar[i])
            mons = monomials(dpol)
            approx_crefs(model)[i] = ConApprox(
                   [integrate(Dict(i => dual(dcon[i]) for i in keys(gmp_variables(model))), 
                              Mom.(mon, jump_function(con))) for mon in mons]
                   .- integrate.(mons, moi_set(con)),
            value.(dvar[i]))
        end    
    end
    return nothing
end

function moments_variable(model::JuMP.Model, v::GMPVariable, deg::Int)
    X = maxdegree_basis(approx_basis(v.v), variables(v.v), deg)
    c = @variable model [1:length(X)] 
    return measure(c, monomials(X))
end

function approximate!(model::GMPModel, ::AbstractPrimalMode)
    degree = approximation_info(model).degree
    # initiate moments 
    pvar = Dict(i => moments_variable(approximation_model(model), v, degree) for (i,v) in gmp_variables(model))
    # add substitutions

    # add measure condition on moments
    just = Dict{Int, Any}()
    scheme_parts = Dict{Int, Vector{SchemePart}}()
    for (i, v) in gmp_variables(model)
        just[i] = []
        scheme_parts[i] = approximation_scheme(model, v)
        for sp in scheme_parts[i]
            cref = primal_scheme_constraint(approximation_model(model), sp, pvar[i])
            push!(just[i], cref)
        end
    end
    # add constraints
    pcon = Dict{Int, Vector{ConstraintRef}}()
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MomentConstraintShape
            cref = @constraint approximation_model(model) integrate(pvar, jump_function(con)) in moi_set(con)            
            pcon[i] = [cref]
        elseif shape(con) isa MeasureConstraintShape
            # add measure constraints
            refmeas = moi_set(con)
            mons = monomials(maxdegree_basis(approx_basis(refmeas), variables(refmeas), approximation_degree(model))) 


            pcon[i] = @constraint approximation_model(model) sum(c.*(MM.expectation.(mons, pvar[index(v)])) for (c, v) in jump_function(con)) .== integrate.(mons, refmeas) 
        end
    end

    # add objective
    @objective approximation_model(model) objective_sense(model) integrate(pvar, objective_function(model))

    optimize!(approximation_model(model))

    for i in keys(gmp_variables(model))
        p = 0
        for (ju, part) in zip(just[i], scheme_parts[i])
            if set_type(part) == PSDCone
                q =  dot(monomials(part.monomials), dual(ju)*monomials(part.monomials))*part.polynomial
            else
                q = dot(dual(ju), monomials(part.monomials)).*part.polynomial
            end
            p = MA.add!(p, q)
        end
        approx_vrefs(model)[i] = VarApprox(MM.Measure(value.(pvar[i].a), pvar[i].x), convert.(Array{Float64}, value.(just[i])), p)
    end
    for (i, con) in gmp_constraints(model)
        if shape(con) isa MeasureConstraintShape
            refmeas = moi_set(con)
            mons = monomials(maxdegree_basis(approx_basis(refmeas), variables(refmeas), approximation_degree(model)))
            #TODO, does minus depend on optimization sense?
            approx_crefs(model)[i] = ConApprox(value.(pcon[i]).-MOI.constant.(JuMP.moi_set.(constraint_object.(pcon[i]))), -dot(dual.(pcon[i]), mons))

        else
            approx_crefs(model)[i] = ConApprox(value.(pcon[i]) .- MOI.constant(moi_set(con)), dual.(pcon[i]))
        end
    end

    return nothing
end

