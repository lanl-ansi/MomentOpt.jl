abstract type AbstractGMPVariable  <: JuMP.AbstractVariable end
include("formulationtypes.jl")
incldue("MeasureJuMP/approximationtypes.jl")
include("MeasureJuMP/measures.jl")
include("MeasureJuMP/continuous.jl")

struct MeasureVariable{T <: AbstractGMPMeasure} <: AbstractGMPVariable
    m::T
end

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, m::AbstractGMPMeasure; extra_kwargs...)
    return MeasureVariable(m)
end

struct ContinuousVariable{T <: AbstractGMPContinuous} <: AbstractGMPVariable
    c::T
end

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, m::AbstractGMPContinuous; extra_kwargs...)
    return ContinuousVariable(m)
end
