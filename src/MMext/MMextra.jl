function MM.measure(μ::GMPVariableRef)
    @assert _is_approximated(μ) "Measure is only available after optimize! has been called." 
    return approx_vrefs(owner_model(μ))[index(μ)].value
end

function MM.moment_matrix(vref::GMPVariableRef, X)
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

function MM.moment_matrix(vref::GMPVariableRef; basis = MonomialBasis)
    X = monomials(maxdegree_basis(basis, variables(vref_object(vref)), Int(floor(approximation_info(owner_model(vref)).degree/2))))
    return moment_matrix(vref, X)
end

Base.broadcastable(moms::MM.Measure) = Ref(moms)

function MM.expectation(p::Number, moms::MM.Measure)
    return p*MM.expectation(_mono_one(variables(moms)), moms)
end
function integrate(p::Number, vref::GMPVariableRef)
    return p*MM.expectation(_mono_one(variables(vref)), measure(vref))
end

function integrate(p::AbstractPolynomialLike, vref::GMPVariableRef)
    return MM.expectation(p, measure(vref))
end

function integrate(model::GMPModel, me::MomentExpr)
    integral = 0
    for (c, m) in me
            integral += integrate(c, m)
        end
    return integral
end

integrate(m::GMPModel, me::AffMomentExpr) = constant(me) + integrate(m, expr(me))

function integrate(ref::Dict, me::MomentExpr)
    integral = 0
    for (c, m) in me 
            integral += MM.expectation(c, ref[index(m)])
    end
    return integral
end

integrate(ref::Dict, me::AffMomentExpr) = constant(me) + integrate(ref, expr(me))

