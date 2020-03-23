module MeasureJuMP

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases

include("relaxtypes.jl")

# MOI extension

using MathOptInterface
const MOI = MathOptInterface
include("attributes.jl")

import LinearAlgebra.dot
using SemialgebraicSets
include("measures.jl")
include("continuous.jl")

include("sets.jl")
# JumP extension

using JuMP

include("variable.jl")

end # module
