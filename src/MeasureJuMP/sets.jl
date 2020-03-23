"""
    EqualToMeasure{T} <: AbstractScalarSet

The set consiting of only a single measure.
"""
struct EqualToMeasure{T <: AbstractGMPMeasure} <: AbstractScalarSet
    measure::T
end

MOI.constant(set::EqualToMeasure) = set.measure

"""
    EqualToContinuous{T} <: AbstractScalarSet

The set consisting of only a single continuous function.
"""
struct EqualToFunction{T <: AbstractGMPConsinuous} <: AbstractScalarSet
    continuous::T
end

MOI.constant(set::EqualToContinuous) = set.continuous

const GMPEqualTo = Union{EqualToMeasure, EqualToContinous}

function Base.:(==)(set1::GMPEqualTo, set2::GMPEqualTo)
    return constant(set1) == constant(set2)
end


for (meas, cont) in zip([],[])
    function MOI.dual_set(set::EqualToMeasure{<:meas})
        return cont(get.(set.measure, [Support(), PolynomialVariables(), RelaxationType(), MomentBasis()])...)
    end
    function MOI.dual_set(set::EqualToContinuous{<:cont})
        return meas(get.(set.measure, [Domain(), PolynomialVariables(), RelaxationType(), MomentBasis()])...)
    end

    MOI.dual_set_type(::EqualToMeasure{<:meas}) = cont
    MOI.dual_set_type(::EqualToContinuous{<:cont}) = meas

end


