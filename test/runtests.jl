using Test

using MomentOpt

using DynamicPolynomials
using SemialgebraicSets
using SumOfSquares
using OrderedCollections
using LinearAlgebra

using SCS
factory =  with_optimizer(SCS.Optimizer, verbose = 1, eps = 1e-16)

include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("model.jl")
include("macros.jl")
include("show.jl")
include("relax.jl")
