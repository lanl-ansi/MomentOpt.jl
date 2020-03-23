"""
    AbstractRelaxationMode

Defines how a model should be relaxed.
"""
abstract type AbstractRelaxationMode end
struct DONT_RELAX <: AbstractRelaxationMode end 
struct PRIMAL_RELAX <: AbstractRelaxationMode end 
struct DUAL_RELAX <: AbstractRelaxationMode end

"""
    AbstractRelaxationType

Defines which relaxation is used to compute a relaxtion for a measure. 
"""
abstract type AbstractRelaxationType end
struct NO_RELAXATION <: AbstractRelaxationType end
struct EXACT_RELAXATION <: AbstractRelaxationType end
abstract type CONIC_FORMULATION <: AbstractRelaxationType end
