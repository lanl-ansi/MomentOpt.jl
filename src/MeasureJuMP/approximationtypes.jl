

"""
    AbstractApproximationType

Defines which formulation is used to relax a measure model or strengthen a continuous model.
"""
abstract type AbstractApproximationType end
struct NO_APPROXIMATION <: AbstractApproximationType end
struct EXACT_APPROXIMATION <: AbstractApproximationType end
abstract type AbstractApproximation <: AbstractApproximationType end


#TODO: Define functions for AbstractApproximation subtypes need to implement
