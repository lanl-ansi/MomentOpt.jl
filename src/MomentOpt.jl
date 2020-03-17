module MomentOpt

import Reexport

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases
using MultivariateMoments
const MM = MultivariateMoments

using SemialgebraicSets
using MathOptInterface
using PolyJuMP
Reexport.@reexport using JuMP
const PJ = PolyJuMP
Reexport.@reexport using SumOfSquares
using OrderedCollections
using LinearAlgebra

const MT = Union{Number, AbstractPolynomialLike}

include("momentsequence.jl")
include("riesz.jl")
include("relaxationtypes.jl")
include("meas.jl")
include("momexpr.jl")
include("momcon.jl")
include("model.jl")
include("relax.jl")
include("postproc.jl")
include("macros.jl")

end# module
