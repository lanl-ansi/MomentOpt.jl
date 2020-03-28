# Measure/Continuous constraints
struct MeasureShape{V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPShape 
    variables::Vector{T}
    basis::T
end

struct ContinuousShape{V <: MP.AbstractVariable, T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPShape 
    variables::Vector{T}
    basis::T
end

JuMP.dual_shape(s::MeasureShape) = ContinuousShape(s.variables, s.basis)
JuMP.dual_shape(s::ContinuousShape) = MeasureShape(s.variables, s.basis)

# Moment/Polynomial constraints
struct MomentShape <: AbstractGMPShape end

struct PolynomialShape <: AbstractGMPShape end

JuMP.dual_shape(::MomentShape) = PolynomialShape()
JuMP.dual_shape(::PolynomialShape) = MomentShape()

# Moment/Monomial substitutions
struct MomentSubstitutionShape{T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPShape 
    basis::T
end

struct IntegralSubstitutionShape{T <: Type{<: MB.AbstractPolynomialBasis}} <: AbstractGMPShape 
    basis::T
end

JuMP.dual_shape(s::MomentSubstitutionShape) = IntegralSubstitutionShape(s.basis)
JuMP.dual_shape(s::IntegralSubstitutionShape) = MomentSubstitutionShape(s.basis)

