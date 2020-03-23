# ModelAttributes

"""
    RelaxationMode(N)
    RelaxationMode()

This attribute defines how the model should be relaxed. By default RelaxationMode() corresponds to no relaxation. The model cannot be approximated in this state.
"""
struct RelaxationMode <: AbstractModelAttribute 
end

"""
    RelaxationStatus()

Attribute for the TerminationStatus of the relaxation of the model.
"""
struct RelaxationStatus <: AbstractModelAttribute end
MOI.is_set_by_optimize(::RelaxationStatus) = true

# VariableAttributes

# Attributes for Measures
"""
    PolynomialVariables()

Attribute for the polynomial variables a measure is acting on. 
"""
struct PolynomialVariables <: AbstractVariableAttribute end

"""
    Support()

Attribute for the support of a measure.
"""
struct Support <: AbstractVariableAttribute end

"""
    RelaxationType()

Attribute for the type of relaxation used for a measure.
"""
struct RelaxationType <: AbstractVariableAttribute end

"""
    MomentBasis()

Attribute for the type of MultivariateBases, used to express the moment of a measure.
"""
struct MomentBasis <: AbstractVariableAttribute end

"""
    MomentFunction()

Attribute for the moment function of a measure.
"""
struct MomentFunction <: AbstractVariableAttribute end

# Attributes for Continuous

"""
    Domain()

Attribute for the domain of a continuous function.
"""
struct Domain <: AbstractVariableAttribute end


