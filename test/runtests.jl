using Test
using MomentOpt
const MO = MomentOpt
using DynamicPolynomials
using SCS
import LinearAlgebra: Symmetric


include("objects.jl")
include("defaultmeasures.jl")
include("variables.jl")
include("affexpr.jl")
include("constraints.jl")

include("gmpmodel.jl")
include("approximate.jl")
include("postproc.jl")
