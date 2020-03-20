export MomentObjective

"""
    MomentObjective

Moment objective.
"""
mutable struct MomObj{PT<:MT}
	sense::MOI.OptimizationSense
	obj::MomExpr{PT}
end

function MomObj(sense::MOI.OptimizationSense, obj::Mom{PT}) where PT<:MT
    return MomObj(sense,convert(MomExpr{PT},obj))
end

function MomObj(sense::MOI.OptimizationSense, pol::PT, meas::Measure) where PT<:MT
    return MomObj(sense,Mom(pol,meas))
end

function MomObj(sense::MOI.OptimizationSense, meas::Measure, pol::PT) where PT<:MT
	return MomObj(sense,Mom(pol,meas))
end

function Base.show(io::IO, f::MomObj)
    if f.sense == MOI.MAX_SENSE
        print(io,"Maximize $(f.obj)")
    else
        @assert f.sense == MOI.MIN_SENSE
        print(io,"Minimize $(f.obj)")
    end
end

function measures(f::MomObj)
	return collect(keys(f.obj.momdict))
end

