export relax!, termination_status

"""
relax 
"""

function relax!(gmp::GMPModel,order::Int,optimizer::OptimizerFactory)
    	# construct dual
	gmp.dual = SOSModel(optimizer)


	if typeof(gmp.objective)==EmptyObjective	
		println("Please define an objective")
		return
	elseif gmp.objective.sense ==:Min
		sossense = MOI.MAX_SENSE
		sign1 = 1
		sign2 =-1
	else
		sossense = MOI.MIN_SENSE
		sign1 =-1
		sign2 = 1
	end

	if length(gmp.constraints) == 0
		println("Please define at least one moment constraint")
		return 
	end

	
	drhs=Dict{Measure, Any}()
	for meas in gmp.measures
		drhs[meas] = 0
	end

	for momcon in gmp.constraints
	if momcon.mode ==:eq
		gmp.dref[momcon] = @variable(gmp.dual)
	elseif momcon.mode ==:leq
		if gmp.objective.sense ==:Min
		gmp.dref[momcon] = @variable(gmp.dual, upper_bound=0)
		else
		gmp.dref[momcon] = @variable(gmp.dual, lower_bound=0)
		end
	elseif momcon.mode ==:geq
		if gmp.objective.sense ==:Min
		gmp.dref[momcon] = @variable(gmp.dual, lower_bound=0) 
		else
		gmp.dref[momcon] =@variable(gmp.dual, upper_bound=0)
		end
	end
	for meas in keys(momcon.lhs.momdict)
		drhs[meas] = drhs[meas] + gmp.dref[momcon]*momcon.lhs.momdict[meas]
	end
	end

	# bookkeeping 
	for meas in gmp.measures	
	gmp.cref[meas] =@constraint(gmp.dual,sign1*get(gmp.objective.obj.momdict,meas,0)+sign2*drhs[meas] in meas.cert, domain = meas.supp, maxdegree = 2*order)

	end

        @objective(gmp.dual,sossense, sum(gmp.dref[momcon]*momcon.rhs for momcon in gmp.constraints))

	optimize!(gmp.dual)

	gmp.dstatus = termination_status(gmp.dual)


end

