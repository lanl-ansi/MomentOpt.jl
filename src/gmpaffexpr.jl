const GMPAffExpr{T <: Number, V <: AbstractGMPVariable} = JuMP.GenericAffExpr{T, V}
include("measureaffexpr.jl")
include("continuousaffexpr.jl")
include("substitutionaffexpr.jl")

