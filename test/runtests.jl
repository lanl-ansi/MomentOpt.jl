using Test
using MomentOpt
const MO = MomentOpt
using DynamicPolynomials

include("objects.jl")
include("variables.jl")
include("affexpr.jl")
include("constraints.jl")

include("gmpmodel.jl")
include("approximate.jl")
include("postproc.jl")

#= 
include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("model.jl")
include("macros.jl")
include("show.jl")
include("relax.jl")
=#
