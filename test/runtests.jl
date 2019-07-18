using MomentOpt
using Test
using DynamicPolynomials
using SemialgebraicSets
using SumOfSquares

using CSDP

@testset "MomentOpt Tests" begin

include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("macros.jl")
include("show.jl")
include("relax.jl")

end
