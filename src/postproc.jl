export moments, atomic, graph, dual_value, christoffel, min_val

function JuMP.objective_value(gmp::GMPModel)
	return objective_value(gmp.dual)
end

function MultivariateMoments.moment_matrix(gmp::GMPModel, measure::Measure)
    return moment_matrix(gmp.cref[measure])
end

function SumOfSquares.moments(gmp::GMPModel, measure::Measure)
	return(SumOfSquares.moments(gmp.cref[measure]))
end

function dual_value(gmp::GMPModel, momcon::JuMP.ConstraintRef)
    return JuMP.value(gmp.dref[momcon.index])
end

function JuMP.value(gmp::GMPModel, mom::Mom{<:Number})
	moms = moments(gmp, mom.meas)
    return dot(moms, constantterm(mom.mon, eltype(monomials(moms))))
end

function JuMP.value(gmp::GMPModel, mom::Mom)
	moms = moments(gmp,mom.meas)
    return dot(moms, mom.mon)
end

function atomic(gmp::GMPModel, measure::Measure; tol=1e-3)
	optmeas = extractatoms(moment_matrix(gmp,measure), tol)
	if typeof(optmeas)== Nothing
		println("Could not detect finite support.")
	else
        optimizers = Dict{Int, Vector{Float64}}()
		for i = 1:length(optmeas.atoms)
		optimizers[i] = optmeas.atoms[i].center
		end
		println("Atomic extraction successful.")
		return optimizers
	end
end

function christoffel(gmp::GMPModel, measure::Measure ; regpar = 1e-8)
    M = moment_matrix(gmp,measure)
    eva,eve = eigen(Matrix(M.Q))
    return sum( (dot(eve[:, i] / √(eva[i] + regpar),M.x))^2 for i = 1:length(eva))
end

function min_val(x::Pair{<:Vector{<:MP.AbstractVariable}, <:Vector{<:Number}},
                 poly::AbstractPolynomialLike)
    p = subs(poly, x)
    if nvariables(p) != 1
        error()
    end
    t = first(variables(p))
    set = algebraicset([differentiate(p, t)])
    Y = [sol[1] for sol in set]
    return Y[argmin([p(t => y) for y in Y])]
end

function min_val_slow(x::Pair{<:Vector{<:MP.AbstractVariable}, <:Vector{<:Number}},
                 Y::StepRangeLen,
                 poly::AbstractPolynomialLike)
    p = subs(poly, x)
    if nvariables(p) != 1
        error()
    end
    t = first(variables(p))
    return Y[argmin([p(t => y) for y in Y])]
end
