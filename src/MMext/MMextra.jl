# extending MultivariateMoments.jl for the use of MomentOpt.jl
Reexport.@reexport using MultivariateMoments
const MM = MultivariateMoments

function MM.measure(μ::GMPVariableRef{AbstractGMPMeasure})
    @assert _is_approximated(μ)
    return approx_vrefs(owner_model(μ))[index(μ)].value
end

function MM.moment_matrix(vref::GMPVariableRef{AbstractGMPMeasure}, X)
    function getmom(i, j)
        x = X[i] * X[j]
        for m in moments(measure(vref))
            if monomial(m) == x
                return moment_value(m)
            end
        end
    end
    return MM.MomentMatrix{Float64}(getmom, X)
end

function MM.moment_matrix(vref::GMPVariableRef)
    X = monomials(variables(vref_object(vref)), 0:Int(floor(approximation_info(owner_model(vref)).degree/2))) # TODO should reflect sparsity 
    return moment_matrix(vref, X)
end
