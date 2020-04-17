abstract type AbstractApproximation end
export ApproximationFunction
"""
    ApproximationFunction

Type used to store a function defined on elements of basis and returning a Number.
"""
struct ApproximationFunction{BT <: MB.AbstractPolynomialBasis} <: AbstractApproximation
    func::Function
    vars::Vector{<:MP.AbstractVariable}
    basis::Type{BT}
end

function (x::ApproximationFunction)(p::AbstractPolynomialLike)
    basis = maxdegree_basis(x.basis, x.vars,maxdegree(p))
    coefs, basis = MB.change_basis(p, basis)
    return dot(coefs, x.func.(basis))
end

# Linear operations on ApproximationFunctions
compatible(f1::ApproximationFunction, f2::ApproximationFunction) = f1.basis == f2.basis && f1.vars == f2.vars

function Base.sum(v::Vector{<:ApproximationFunction})
    f1 = first(v)
    for f2 in v[2:end]
        @assert compatible(f1, f2) "Only functions acting on same basis can be added."
    end
    return ApproximationFunction(x -> sum(f.func(x) for f in v), f1.vars, f1.basis)
end

function Base.:+(f1::ApproximationFunction, f2::ApproximationFunction)
    return sum([f1, f2])
end
Base.:*(a::Number, f::ApproximationFunction) = ApproximationFunction(x -> a*f.func(x), f.vars, f.basis)
Base.:*(f::ApproximationFunction, a::Number) = a*f
Base.:/(f::ApproximationFunction, a::Number) = inv(a)*f
Base.:-(f::ApproximationFunction) = (-1)*f
Base.:-(f1::ApproximationFunction, f2::ApproximationFunction) = f1 + (-f2)


"""
    AbstractApproximationScheme

Abstract type definines how a GMPObject should be approximated. Every implementation T of AbstractApproximationScheme should define a function approximation_scheme(scheme::T, K::AbstractSemialgebraicSet, d::Int) returning a dictionary explaining the scheme. If sparsity or symmetry can be taken into account the function sparse_approximation_scheme(scheme::T, f::AbstractPolynomialLike, K::AbstractSemialgebraicSet, d::Int) invariant_approximation_scheme(scheme::T, action, K::AbstractSemialgebraicSet, d::Int) should be defined. Otherwise those functions defautlt to approximation_scheme.
"""
abstract type AbstractApproximationScheme <: AbstractApproximation end

function approximation_scheme(scheme::AbstractApproximationScheme, K::AbstractSemialgebraicSet, d::Int) end

sparse_approximation_scheme(scheme::AbstractApproximationScheme, f::MP.AbstractPolynomialLike, K::AbstractSemialgebraicSet, d::Int) = approximation_scheme(scheme, K, d)

sparse_approximation_scheme(scheme::AbstractApproximationScheme, action, K::AbstractSemialgebraicSet, d::Int) = approximation_scheme(scheme, K, d)


