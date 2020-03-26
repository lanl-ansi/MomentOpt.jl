"""
    EqualToMeasure{T} <: AbstractScalarSet

The set consiting of only a single measure.
"""
struct EqualToMeasure{T <: AbstractGMPMeasure} <: MOI.AbstractScalarSet
    measure::T
end

MOI.constant(set::EqualToMeasure) = set.measure
MOI.dual_set(s::EqualToMeasure) = GreaterThanContinuous(VariableContinuous(get.(constant(s), GenericMeasureAttributes))...)
MOI.dual_set_type(::Type{EqualToMeasure}) = GreaterThanContinuous

"""
    GreaterThanContinuous{T} <: AbstractScalarSet

The set consisting of all continuous functions bigger than a single continuous.
"""
struct GreaterThanContinuous{T <: AbstractGMPContinuous} <: MOI.AbstractScalarSet
    continuous::T
end

MOI.constant(set::GreaterThanContinuous) = set.continuous
MOI.dual_set(s::GreaterThanContinuous) = EqualToMeasure(VariableMeasure(get.(constant(s), GenericContinuousAttributes))...)
MOI.dual_set_type(::Type{GreaterThanContinuous}) = EqualToMeasure

const GMPset = Union{EqualToMeasure, GreaterThanContinuous}

function Base.:(==)(set1::GMPset, set2::GMPset)
    return constant(set1) == constant(set2)
end
