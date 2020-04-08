@testset "AbstractGMPExpr" begin
# gmp_coefficients(e::AbstractGMPExpr) 
# gmp_variables(e::AbstractGMPExpr)
# Base.iszero(e::AbstractGMPExpr)
# Base.getindex(e::AbstractGMPExpr, s) 
# Base.length(e::AbstractGMPExpr) 
# Base.iterate(e::AbstractGMPExpr)
# Base.iterate(e::AbstractGMPExpr, s)
# Base.isequal(me1::AbstractGMPExpr, me2::AbstractGMPExpr)
@test MO.GMPEmptyExpr() isa MO.AbstractGMPExpressionLike
end

@testset "AbstractGMPExpr - Arith" begin
    @testset "ObjectExpr - Arith" begin

    end
    @testset "MomentExpr - Arith" begin

    end

end

@testset "GMPAffExpr - Arith" begin

    @testset "AffObjectExpr - Arith" begin

    end

end

@testset "ObjectExpr" begin
 # same_variables(vrefs::Vector{GMPVariableRef})
 # same_vref_type(vrefs::Vector{GMPVariableRef})
# ObjectExpr(c::T, v::GMPVariableRef) where {T}
# JuMP.function_string(io, e::ObjectExpr)
# MP.variables(e::ObjectExpr) 
# vref_type(e::ObjectExpr) 
# degree(e::ObjectExpr) 
end


@testset "AffObjectExpr" begin

end

@testset "MomentExpr" begin
# MomentExpr(cv::Vector{<:Union{Number, MP.AbstractPolynomialLike}}, mv::Vector{GMPVariableRef})
# MomentExpr(cv::Union{Number, MP.AbstractPolynomialLike}, mv::Union{GMPVariableRef, ObjectExpr})# degree(e::MomentExpr) = maximum(MP.maxdegree.(gmp_coefficients(e)))
# function JuMP.function_string(io, e::MomentExpr)
# momexp_by_measure(e::MomentExpr{T, S}) 

end

