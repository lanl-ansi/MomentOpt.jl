@testset "Approximate" begin

    @testset "ApproximationSequence" begin
        # MO.gmp_object(ms::ApproximationSequence)
        # Base.show(io::IO, ms::ApproximationSequence)
        # Base.length(ms::ApproximationSequence)
        # Base.getindex(ms::ApproximationSequence, m::MP.AbstractPolynomialLike)
        # Base.getindex(ms::ApproximationSequence, m::Number)
        # Base.getindex(ms::ApproximationSequence, v::AbstractArray)
        # MO.approximation_sequence(model::JuMP.Model, o::AbstractGMPObject, degree::Int)
    end
    @testset "Riesz" begin
        # Base.show(io::IO, ap::GMPObjectApproximation)
        # approximation_sequence(model::GMPModel, index::Int)
        # riesz(m::GMPModel, index::Int, p::MP.AbstractPolynomialLike)
        # riesz(ms::GMPModel, index::Int, m::Number)
        # riesz(model::GMPModel, me::MomentExpr)
        # riesz(ms::GMPModel, me::AffMomentExpr) =
    end

    @testset "StandardRelaxation - Primal" begin
        # MO.QuadraticModuleMonomials(o::AbstractGMPObject, degree::Int)
        # MO.MomentPutinar(model::GMPModel, degree::Int)

    end

    @testset "StandardRelaxation - Dual" begin

    end

end
