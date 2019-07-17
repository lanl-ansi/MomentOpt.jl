module MomentOpt 

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets
using MathOptInterface
using PolyJuMP
using JuMP
const PJ = PolyJuMP
using SumOfSquares

const MT = Union{Number,AbstractPolynomialLike}

include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("model.jl")
include("relax.jl")
include("postproc.jl")
include("macros.jl")


end# module
