export moments, atomic, graph, dual_value

function JuMP.objective_value(gmp::GMPModel)
	return objective_value(gmp.dual)
end

function SumOfSquares.moments(gmp::GMPModel, measure::Measure)
	return(SumOfSquares.moments(gmp.cref[measure]))
end

function dual_value(gmp::GMPModel,momcon::MomCon)
	return JuMP.value(gmp.dref[momcon])
end

function dual_value(gmp::GMPModel,momcons::Vector{MomCon})
	dv = Float64[]
	for momcon in momcons
		push!(dv,JuMP.value(gmp.dref[momcon]))
	end
	return dv
end

function JuMP.value(gmp::GMPModel, mom::Mom{T}) where T<:Number
	moms = moments(gmp,mom.meas)
	idx = findfirst(x->x==1,moms.x)
	return mom.mon*moms.a[idx]
end

function JuMP.value(gmp::GMPModel,mom::Mom{T}) where T<:AbstractMonomialLike
	moms = moments(gmp,mom.meas)
	idx = findfirst(x->x==mom.mon,moms.x)
	return moms.a[idx]
end

function JuMP.value(gmp::GMPModel,mom::Mom{T}) where T<:AbstractPolynomialLike
	moms = moments(gmp,mom.meas)
	vec = []
	for mon in mom.mon.x
		idx=findfirst(x->x==mon,moms.x)
		push!(vec,moms.a[idx])
	end
	return sum(mom.mon.a[i]*vec[i] for i = 1:length(mom.mon.a))
end

function JuMP.value(gmp::GMPModel, vec::Vector{T}) where T<:AbstractMoment
	val = []
	for mom in vec
	push!(val,value(gmp,mom))
	end
	return val
end


function atomic(gmp::GMPModel, measure::Measure)
	optmeas = extractatoms(moment_matrix(gmp.cref[measure]),1e-03)
	if typeof(optmeas)== Nothing
		println("Could not detect finite support.")
	else
		optimizers = Dict()
		for i = 1:length(optmeas.atoms)
		optimizers[i] = optmeas.atoms[i].center
		end
		println("Atomic extraction successful.")
		return optimizers
	end
end

function graph(gmp::GMPModel, measure::Measure)
println("Coming soon...")
end


