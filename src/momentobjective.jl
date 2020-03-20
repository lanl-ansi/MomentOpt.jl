export MomentObjective

"""
    MomentObjective

Type holding a moment objective.
"""
mutable struct MomentObjective{T <: Number, PT <: MT}
	sense::MOI.OptimizationSense
	obj::MomentExpr{T, PT}
end

JuMP.objective_function(mo::MomentObjective) = mo.obj
JuMP.objective_sense(mo::MomentObjective) = mo.sense
measures(f::MomentObjective) = measures(objective_function(mo))

function Base.show(io::IO, f::MomentObjective)
    if f.sense == MOI.MAX_SENSE
        print(io,"Maximize $(f.obj)")
    elseif f.sense == MOI.MIN_SENSE
        print(io,"Minimize $(f.obj)")
    else
        nothing
    end
end
