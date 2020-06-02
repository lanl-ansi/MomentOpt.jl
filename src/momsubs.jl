"""
    MomentSubstitute

"""
struct MomentSubstitute{T} <: MOI.AbstractScalarSet
    momexpr::T
end

MOI.constant(ms::MomentSubstitute) = ms.momexpr

function Base.:(==)(set1::MomentSubstitute{T}, set2::MomentSubstitute{T}) where T 
    return constant(set1) == constant(set2)
end
JuMP.sense_to_set(::Function, ::Val{:(>>)}) = MomentSubstitute(zero(MomentExpr{Int, Int}))
MOIU.shift_constant(set::MomentSubstitute, value) = MomentSubstitute(MOI.constant(set) + value)


"""
    MomentSubstitutionShape 

"""
struct MomentSubstitutionShape <: AbstractGMPShape end

JuMP.reshape_vector(expr::MomentExpr, ::MomentSubstitutionShape) = expr
JuMP.reshape_set(set::MOI.AbstractScalarSet, ::MomentSubstitutionShape) = set

"""
    MomentSubstitution{C, T} <: AbstractGMPConstraint

"""
struct MomentSubstitution{C, T} <: AbstractGMPConstraint
    func::C
    set::MomentSubstitute{T}
end

JuMP.jump_function(con::MomentSubstitution) = con.func
JuMP.moi_set(con::MomentSubstitution) = con.set
JuMP.shape(con::MomentSubstitution) = MomentSubstitutionShape()


function JuMP.build_constraint(_error::Function, func::MomentExpr,
                               set::MomentSubstitute)

    # check that func is an actual moment
    # check that set does not involve a moment in the measure of func with higher degree
        
    return MomentSubstitution(func, set)
end
