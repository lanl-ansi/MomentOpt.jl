using Test

using MomentOpt

using DynamicPolynomials
using SemialgebraicSets
using SumOfSquares
using OrderedCollections
using LinearAlgebra

using SCS
factory = SCS.Optimizer

include("momentsequence.jl")
include("riesz.jl")
include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("model.jl")
include("macros.jl")
include("show.jl")
include("relax.jl")
