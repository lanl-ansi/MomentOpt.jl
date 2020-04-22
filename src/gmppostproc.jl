function JuMP.objective_value(model::GMPModel) 
    if approximation_mode(model) isa AbstractPrimalMode
        return objective_value(approximation_model(model))
    elseif  approximation_mode(model) isa AbstractDualMode
        return dual_objective_value(approximation_model(model))
    end
end
function JuMP.dual_objective_value(model::GMPModel) 
    if approximation_mode(model) isa AbstractPrimalMode
        return dual_objective_value(approximation_model(model))
    elseif  approximation_mode(model) isa AbstractDualMode
        return objective_value(approximation_model(model))
    end
end
function JuMP.primal_status(model::GMPModel) 
    if approximation_mode(model) isa AbstractPrimalMode
        return primal_status(approximation_model(model))
    elseif  approximation_mode(model) isa AbstractDualMode
        return dual_status(approximation_model(model))
    end
end
function JuMP.dual_status(model::GMPModel) 
    if approximation_mode(model) isa AbstractPrimalMode
        return dual_status(approximation_model(model))
    elseif  approximation_mode(model) isa AbstractDualMode
        return primal_status(approximation_model(model))
    end
end

function JuMP.termination_status(model::GMPModel)
    return termination_status(approximation_model(model))
end

export atomic
"""
    atomic(vref::GMPVariableRef; tol=1e-3)

Tries to exctract the atomic support of a measure. 
"""
function atomic(vref::GMPVariableRef; tol = 1e-3)
	optmeas = extractatoms(moment_matrix(vref), tol)
	if typeof(optmeas) == Nothing
        return nothing
	else
        optimizers = Dict{Int, Vector{Float64}}()
		for i = 1:length(optmeas.atoms)
		    optimizers[i] = optmeas.atoms[i].center
		end
		return optimizers
	end
end


function min_val(x::Pair{<:Vector{<:MP.AbstractVariable}, <:Vector{<:Number}},
                 poly::AbstractPolynomialLike)
    p = subs(poly, x)
    t = first(variables(p))
    set = algebraicset([differentiate(p, t)])
    Y = [sol[1] for sol in set]
    return Y[argmin([p(t => y) for y in Y])]
end

function christoffel(vref::GMPVariableRef; regpar = 1e-8)
    M = moment_matrix(vref)
    eva, eve = eigen(MM.getmat(M))
    return sum( (dot(eve[:, i] / √(eva[i] + regpar), M.basis.monomials))^2 for i = 1:length(eva))
end

export graph
"""
    graph(vref::GMPVariable, indep_vars::Vector{MP.AbstractVariable})

    Tries to extract the graph of a scalar function f such that [indep_vars, f(indep_vars)] is the support of vref.
"""
function graph(vref::GMPVariableRef, indep_vars::Vector{<:MP.AbstractVariable}; regpar = 1e-8)
    ch =  christoffel(vref; regpar = regpar)
    @assert length(setdiff(variables(ch), indep_vars)) == 1 "Only scalar functions supported for graph."
    return x -> min_val(indep_vars => x, ch)
end
