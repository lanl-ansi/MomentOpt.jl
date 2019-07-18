using MomentOpt
using Test
using DynamicPolynomials
using SemialgebraicSets
using SumOfSquares


@testset "MomentOpt Tests" begin

include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("macros.jl")
include("show.jl")

end
