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


Base.broadcastable(m::MM.Measure) = Ref(m)

function integrate(p::Number, moms::MM.Measure)
    @assert moms.x[end] == 1
    return p * moms.a[end]
end

function integrate(p::MP.AbstractPolynomialLike, moms::MM.Measure)
    integral = 0.0
    for (c, m) in zip(coefficients(p), monomials(p))
        idx = findfirst(x -> x==m, monomials(moms))
        integral += c*moms.a[idx]
    end
    return integral
end

function integrate(model::GMPModel, me::MomentExpr)
    integral = 0
    for (c, m) in momexp_by_measure(me)
        for (k, meas) in m
            integral += integrate(c*k, approx_vrefs(model)[index(meas)].value)
        end
    end
    return integral
end

integrate(m::GMPModel, me::AffMomentExpr) = constant(me) + integrate(m, expr(me))
