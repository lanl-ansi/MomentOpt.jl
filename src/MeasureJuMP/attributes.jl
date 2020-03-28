# ModelAttributes

"""
    ApproximationMode()

Attribute for the mode used to approximate a model.
"""
struct ApproximationMode <: MOI.AbstractModelAttribute end

"""
    ApproximationStatus()

Attribute for the TerminationStatus of the approximation of the model.
"""
struct ApproximationStatus <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::ApproximationStatus) = true
