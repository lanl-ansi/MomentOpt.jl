export moments, atomic, graph, dual_value

function JuMP.objective_value(gmp::GMPModel)
	return objective_value(gmp.dual)
end

function moments(gmp::GMPModel, measure::Measure)
	return(moment_matrix(gmp.cref[measure]))
	# should be replaced by the actual command
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


