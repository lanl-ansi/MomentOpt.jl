# extending MultivariateMoments.jl for the use of MomentOpt.jl
Reexport.@reexport using MultivariateMoments
const MM = MultivariateMoments

struct GMPMoment{T} <: MM.AbstractMoment{T}
    α::T
    x::MP.AbstractPolynomialLike
    meas::GMPVariableRef
end

"""
    moment(monom::AbstractPolynomialLike, meas::GMPVariableRef; value = )

Creates a moment from the input data.
"""
function moment(monom::MP.AbstractPolynomialLike, meas::GMPVariableRef; value = approximation(meas, monom))
    @assert vref_object(meas) isa AbstractGMPMeasure "GMPVariableRef has to refer to an AbstractGMPMeasure."
    return GMPMoment(value, monom, meas)
end

MM.moment_value(m::GMPMoment) = m.α
MP.monomial(m::GMPMoment) = m.x

function MM.moments(μ::GMPVariableRef) 
    dict = approximation(μ)
    return map((k, v) -> moment(k, μ; value = v), keys(dict), values(dict))
end

function MM.moment_matrix(vref::GMPVariableRef, X)
    function getmom(i, j)
        x = X[i] * X[j]
        for m in moments(vref)
            if monomial(m) == x
                return moment_value(m)
            end
        end
        throw(ArgumentError("Measure does not have the moment $(x)"))
    end
    return MM.MomentMatrix{Float64}(getmom, X)
end

function MM.moment_matrix(vref::GMPVariableRef)
    X = monomials(MOI.get(vref_object(vref), Variables()), 0:Int(floor(approximation_info(owner_model(vref)).degree/2))) # TODO should reflect sparsity 
    return moment_matrix(vref, X)
end

