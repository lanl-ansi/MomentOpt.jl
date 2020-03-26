# ModelAttributes

"""
    ApproximationMode()

Attribute for the mode used to approximate a model.
"""
struct ApproximationMode <: AbstractModelAttribute 
end

"""
    ApproximationStatus()

Attribute for the TerminationStatus of the approximation of the model.
"""
struct ApproximationStatus <: AbstractModelAttribute end
MOI.is_set_by_optimize(::ApproximationStatus) = true

# VariableAttributes

"""
    Variables()

Attribute for the variables a measure is acting on. 
"""
struct Variables <: AbstractVariableAttribute end

"""
    Support()

Attribute for the support of a measure.
"""
struct Support <: AbstractVariableAttribute end

"""
    ApproximationType()

Attribute for the type of approximation used for a measure/continuous function.
"""
struct ApproximationType <: AbstractVariableAttribute end

"""
    GMPBasis()

Attribute for the type of MultivariateBases, used to express the moment of a measure or the monomials of a polynomial.
"""
struct GMPBasis <: AbstractVariableAttribute end

"""
    MomentFunction()

Attribute for the moment function of a measure.
"""
struct MomentFunction <: AbstractVariableAttribute end

"""
    Domain()

Attribute for the domain of a continuous function.
"""
struct Domain <: AbstractVariableAttribute end
"""
    CoefFunction()

Attribute for the coefficient function of a continuous function.
"""
struct CoefFunction <: AbstractVariableAttribute end

const GenericMeasureAttributes = [Variables(), Support(), ApproximationType(), GMPBasis()]
const GenericContinuousAttributes = [Variables(), Domain(), ApproximationType(), GMPBasis()] 
