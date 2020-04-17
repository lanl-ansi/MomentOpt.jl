abstract type AbstractGMPObject end
Base.broadcastable(o::AbstractGMPObject) = Ref(o)

# AbstractGMPObject functions

MP.variables(m::AbstractGMPObject) = m.variables
approx_basis(m::AbstractGMPObject) = m.approx_basis

"""
    covering_basis(t::AbstractGMPObject, p::MP.AbstractPolynomialLike)
    covering_basis(t::AbstractGMPObject, p::Number)

Returns all monomials in the basis of t covering the monomials of p.    
"""
function covering_basis(t::AbstractGMPObject, p::MP.AbstractPolynomialLike)
    @assert variables(p) ⊆  variables(t) "Object does not act on $(variables(p))."
    return MB.basis_covering_monomials(approx_basis(t), monomials(p)) 
end

"""
    _mono_one(vars)

returns the monomial 1 defined on vars. 
"""
_mono_one(vars::Vector{<:MP.AbstractVariable}) = prod(var^0 for var in vars)

covering_basis(t::AbstractGMPObject, p::Union{Number, AbstractJuMPScalar}) = covering_basis(t, p*_mono_one(variables(t)))

"""
    maxdegree_basis(t::AbstractGMPObject, d::Int)

Returns all monomials up to degree d in the basis of t.    
"""
function MB.maxdegree_basis(t::AbstractGMPObject, d::Int)
    return maxdegree_basis(approx_basis(t), variables(t), d) 
end

"""
    AbstractGMPMeasure

Abstract type to represent measures in MomentOpt. Implementation of subtypes should implement MOI.get functions for `::Variables`, `::Support`, `::ApproximationType`, and `::GMPBasis`.
if get(m, ApproximationType()) == EXACT_APPROXIMATION, the measure m should have a field m.moment_function::Function. 
By default there are four implementations of AbstractGMPMeasure:
  * [AnalyticMeasure](@ref) :
  * [VariableMeasure](@ref) :

"""
abstract type AbstractGMPMeasure <: AbstractGMPObject end

export AnalyticMeasure

"""
    AnalyticMeasure{V, T} <: AbstractGMPMeasure

Type representing an analytic measure. Its field approx_function allows to integrate polynomials."""
struct AnalyticMeasure{V <: MP.AbstractVariable, BT <: MB.AbstractPolynomialBasis} <: AbstractGMPMeasure
    variables::Vector{V}
    approx_basis::Type{BT}
    approx_function::ApproximationFunction{BT}
end

function AnalyticMeasure(variables::Vector{V}, moment_basis::T, moment_function::Function) where  {V <: MP.AbstractVariable, T <: Type{<:MB.AbstractPolynomialBasis}} 
    return AnalyticMeasure(variables, moment_basis, ApproximationFunction(moment_function, variables, moment_basis))
end

constructor(::Type{AnalyticMeasure{S,T}}) where {S, T} = (x, y, z) -> AnalyticMeasure(x, y, z)
Base.show(io::IO, ::AnalyticMeasure) = print(io, "AnalyticMeasure")

export VariableMeasure

"""
    VariableMeasure{V, S, R, T}

Type representing a measure that can be relaxed via a conic relaxation. 
"""
struct VariableMeasure{V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet, BT <: MB.AbstractPolynomialBasis, R <: AbstractApproximationScheme} <: AbstractGMPMeasure
    variables::Vector{V}
    bsa_set::S
    approx_basis::Type{BT}
    approx_scheme::R
end

export Meas
"""
    Meas(vars; support = FullSpace(), basis = MonomialBasis, scheme = PutinarScheme())

Shortcut to define a VariableMeasure.
"""
Meas(vars; support = FullSpace(), basis = MonomialBasis, scheme = PutinarScheme()) = VariableMeasure(vars, support, basis, scheme)

"""
    AbstractGMPContinuous

Under some regularity assumptions one can say that the dual space of measures is the space of continuous functions. This is for example the case if the support of the measure is compact. This abstract type allows to define custom testfunctions.
"""
abstract type AbstractGMPContinuous <: AbstractGMPObject end

export AnalyticContinuous

"""
    AnalyticContinuous{V, T} <: AbstractGMPContinuous

Type representing an analytic continuous function. Its field coef_function allows to compute the coefficients for arbitrary close polynomial approximations.
"""
struct AnalyticContinuous{V <: MP.AbstractVariable, BT <: MB.AbstractPolynomialBasis} <: AbstractGMPContinuous
    variables::Vector{V}
    approx_basis::Type{BT}
    approx_function::ApproximationFunction{BT}
end

function AnalyticContinuous(variables::Vector{V}, monom_basis, monom_function::Function) where {V <: MP.AbstractVariable} 
    return AnalyticContinuous(variables, monom_basis, ApproximationFunction(monom_function, variables, monom_basis)) 
end

Base.show(io::IO, ::AnalyticContinuous) = print(io, "AnalyticContinuous")
constructor(::Type{AnalyticContinuous{S,T}}) where {S, T} = (x, y, z) -> AnayticContinuous(x, y, z)

export VariableContinuous

"""
    VariableContinuous{V, S, R}

Type representing a continuous funciton that can be stregthend via a conic approximation. 
"""
struct VariableContinuous{V <: MP.AbstractVariable, S <: AbstractBasicSemialgebraicSet, BT <: MB.AbstractPolynomialBasis, R <: AbstractApproximationScheme} <: AbstractGMPContinuous
    variables::Vector{V}
    bsa_set::S
    approx_basis::Type{BT}
    approx_scheme::R
end

export Cont
"""
    Meas(vars; support = FullSpace(), basis = MonomialBasis, scheme = PutinarScheme())

Shortcut to define a VariableMeasure.
"""
Cont(vars; domain = FullSpace(), basis = MonomialBasis, scheme = PutinarScheme()) = VariableContinuous(vars, domain, basis, scheme)

#
const AnalyticGMPObject = Union{AnalyticMeasure, AnalyticContinuous}
# we can perform linear operations with analytic measures by performing them on the respective approximation funciton
#
compatible(f1::AnalyticGMPObject, f2::AnalyticGMPObject) = variables(f1) == variables(f2) && supertype(typeof(f1)) == supertype(typeof(f2))
function Base.sum(v::Vector{<:AnalyticGMPObject})
    f1 = first(v)
    for f2 in v[2:end]
        @assert compatible(f1, f2) "Cannot add GMPObjects, because they are either acting on different variables or are of different Type."
    end
    return constructor(typeof(v[1]))(variables(v[1]), approx_basis(v[1]), sum(approx_function.(v))) 
end
function Base.:+(f1::AnalyticGMPObject, f2::AnalyticGMPObject)
@assert compatible(f1, f2) "Cannot add GMPObjects, because they are either acting on different variables or are of different Type."
    return sum([f1, f2])
end
Base.:*(a::Number, f::AnalyticGMPObject) = constructor(typeof(f))(variables(f), approx_basis(f), a*approx_function(f)) 
Base.:*(f::AnalyticGMPObject, a::Number) = a*f
Base.:/(f::AnalyticGMPObject, a::Number) = inv(a)*f
Base.:-(f::AnalyticGMPObject) = (-1)*f
Base.:-(f1::AnalyticGMPObject, f2::AnalyticGMPObject) = f1 + (-f2)

# functions only available for AnalyticGMPObjects

approx_function(m::AnalyticGMPObject) = m.approx_function

export integrate, approx

function approx(p::Union{Number, AbstractJuMPScalar}, m::AnalyticGMPObject)
    return p*approx_function(m)(first(monomials(covering_basis(m, p))))
end

function approx(p::MP.AbstractPolynomialLike, m::AbstractGMPMeasure)
    @assert variables(p) ⊆ variables(m)
    return approx_function(m)(p)
end

"""
     integrate(p::Number, m::AnalyticMeasure)
     integrate(p::MP.AbstractPolynomialLike, m::AnalyticMeasure)

Returns the integral of p with respect to m. 
"""
const integrate = approx

export partial_integrate, partial_approx

function partial_approx(p::MP.AbstractPolynomialLike, m::AnalyticGMPObject)
    if variables(p) ⊆  variables(m)
        return approx(p, m)
    else
       coefs, mons = polypoly(p, variables(m))
       return dot(coefs, approx.(mons, m))
    end
end


"""
     partial_integrate(p::MP.AbstractPolynomialLike, m::AbstractGMPMeasure)

Returns the integral of p with respect to m, in case m does not cover all variables of p. The result is a polynomial in the remaining variables.  
"""
const partial_integrate = partial_approx

#
const VariableGMPObject = Union{VariableMeasure, VariableContinuous}
# functions only available for VariableGMPObjects

bsa_set(v::VariableGMPObject) = v.bsa_set
approx_scheme(v::VariableGMPObject) = v.approx_scheme
