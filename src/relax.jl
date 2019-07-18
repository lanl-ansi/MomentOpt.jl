export relax!, termination_status

function add_prod!(p::AbstractPolynomial{JuMP.AffExpr}, q::AbstractPolynomialLike{<:Number}, args...)
    pt = terms(p)
    qt = terms(q)
    pelst = iterate(pt)
    qelst = iterate(qt)
    r = zero(q)
    while qelst !== nothing
        qel, qst = qelst
        if pelst === nothing
            r += qel
            qelst = iterate(qt, qst)
        else
            pel, pst = pelst
            if monomial(pel) == monomial(qel)
                JuMP.add_to_expression!(coefficient(pel), coefficient(qel), args...)
                qelst = iterate(qt, qst)
                pelst = iterate(pt, pst)
            elseif monomial(pel) < monomial(qel)
                r += qel
                qelst = iterate(qt, qst)
            else
                pelst = iterate(pt, pst)
            end
        end
    end
    if iszero(nterms(r))
        return p
    else
        return p + *(r, args...)
    end
end
function add_prod!(p::AbstractPolynomial{JuMP.AffExpr}, q::Number, args...)
    add_prod!(p, constantterm(q, p), args...)
end



"""
    relax!(gmp::GMPModel, order::Int, optimizer::OptimizerFactory)

Computes the a relaxation of the generlized moment problem according to the certificates specified for each measure. 
"""
function relax!(gmp::GMPModel, order::Int, optimizer::OptimizerFactory)
    	# construct dual
	gmp.dual = SOSModel(optimizer)


	if typeof(gmp.objective) == EmptyObjective	
		println("Please define an objective")
		return
	elseif gmp.objective.sense == MOI.MIN_SENSE
		sossense = MOI.MAX_SENSE
		sign1 =  1
	else
		sossense = MOI.MIN_SENSE
		sign1 = -1
	end

	if isempty(gmp.constraints)
		@error "Define at least one moment constraint!"
	end
	
    PT = polynomialtype(variabletype(first(gmp.measures)), JuMP.AffExpr)
    drhs = Dict{Measure, PT}()
	for meas in gmp.measures
        drhs[meas] = zero(PT)
	end

    for momcon in gmp.constraints
        if momcon.set isa MOI.EqualTo
            gmp.dref[momcon] = @variable(gmp.dual)
        elseif momcon.set isa MOI.LessThan
            if gmp.objective.sense == MOI.MIN_SENSE
                gmp.dref[momcon] = @variable(gmp.dual, upper_bound=0)
            else
                gmp.dref[momcon] = @variable(gmp.dual, lower_bound=0)
            end
        elseif momcon.set isa MOI.GreaterThan
            if gmp.objective.sense == MOI.MIN_SENSE
                gmp.dref[momcon] = @variable(gmp.dual, lower_bound=0) 
            else
                gmp.dref[momcon] = @variable(gmp.dual, upper_bound=0)
            end
        end
        for meas in measures(momcon)
            drhs[meas] = add_prod!(drhs[meas], momcon.func.momdict[meas], gmp.dref[momcon])
        end
    end

    # bookkeeping 
	for meas in gmp.measures
        p = drhs[meas]
        if sossense == MOI.MAX_SENSE
            p = -p
        end
        obj = get(gmp.objective.obj.momdict ,meas, nothing)
        if obj != nothing
            if sossense == MOI.MIN_SENSE
                obj = -obj
            end
            # Modifies drhs[meas] but it's ok since it not used anywhere else
            p = add_prod!(p, obj)
        end
	    gmp.cref[meas] = @constraint(gmp.dual, p in meas.cert, domain = meas.supp, maxdegree = 2*order)
	end

        @objective(gmp.dual,sossense, sum(gmp.dref[momcon]*momcon.rhs for momcon in gmp.constraints))

	optimize!(gmp.dual)

	gmp.dstatus = termination_status(gmp.dual)
end
