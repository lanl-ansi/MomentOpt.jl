export moments, atomic, graph, dual_value, christoffel, min_val
function MultivariateMoments.moment_matrix(gmp::GMPModel, vref::GMPVariableRef{AbstractGMPMeasure})
    @assert is_approximated(vref)
    return moment_matrix(approx_vrefs[index(vref)])
end

function atomic(gmp::GMPModel, measure::Measure; tol=1e-3, print_level = 1)
	optmeas = extractatoms(moment_matrix(gmp,measure), tol)
	if typeof(optmeas)== Nothing
        if print_level ==1
		    println("Could not detect finite support.")
        end
	else
        optimizers = Dict{Int, Vector{Float64}}()
		for i = 1:length(optmeas.atoms)
		optimizers[i] = optmeas.atoms[i].center
		end
        if print_level ==1 
	    	println("Atomic extraction successful.")
        end
		return optimizers
	end
end

function christoffel(gmp::GMPModel, measure::Measure ; regpar = 1e-8)
    M = moment_matrix(gmp,measure)
    eva,eve = eigen(Matrix(M.Q))
    return sum( (dot(eve[:, i] / √(eva[i] + regpar), M.basis.monomials))^2 for i = 1:length(eva))
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
