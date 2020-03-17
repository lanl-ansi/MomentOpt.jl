

abstract type AbstractConicFormulation end
abstract type AbstractPrimalFormulation <: AbstractConicFormulation end
abstract type AbstractDualFormulation <: AbstractConicFormulation end


struct PrimalPutinar <: AbstractPrimalFormulation 
end

struct DualPutinar <: AbstractDualFormulation
end
