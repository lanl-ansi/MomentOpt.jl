abstract type AbstractApproximationScheme end


export KeepEquality, ReduceEquality, SplitEquality
abstract type AbstractEqualityHandler end
struct KeepEquality <: AbstractEqualityHandler end
struct ReduceEquality <: AbstractEqualityHandler end
struct SplitEquality <: AbstractEqualityHandler end

struct SchemePart{T<:Union{JuMP.PSDCone, MOI.AbstractSet}}
    polynomial::MP.AbstractPolynomialLike
    monomials::MB.AbstractPolynomialBasis
    moi_set::T
end

set_type(::SchemePart{T}) where T = T

approximation_scheme(model::JuMP.AbstractModel, v) = approximation_scheme(approx_scheme(v.v), bsa_set(v.v), variables(v.v), approximation_degree(model))

function dual_scheme_variable(model::JuMP.Model, sp::SchemePart{<:MOI.Nonnegatives})
    vref = @variable model lower_bound = 0
    return vref
end

function dual_scheme_variable(model::JuMP.Model, sp::SchemePart{<:MOI.Zeros})
    vref = @variable model variable_type = Poly(sp.monomials)
    return vref
end

function dual_scheme_variable(model::JuMP.Model, sp::SchemePart{PSDCone})
    vref = @variable model variable_type = SOSPoly(sp.monomials)
    return vref
end

function primal_scheme_part(sp::SchemePart{<:Union{MOI.Nonnegatives, MOI.Zeros}}, moment_vars::MM.Measure)
    return MM.expectation.(monomials(sp.monomials)*sp.polynomial, moment_vars)
end

function localization_matrix(mons::MB.AbstractPolynomialBasis, poly::MP.AbstractPolynomialLike, moment_vars::MM.Measure)
    data = zeros(JuMP.AffExpr, length(mons), length(mons))
    mons = monomials(mons)
    for (i, moni) in enumerate(mons)
        for j in i:length(mons)
            data[i, j] = MM.expectation(moni*mons[j]*poly, moment_vars)
        end
    end
    return Symmetric(data)
end

function primal_scheme_part(sp::SchemePart{PSDCone}, moment_vars::MM.Measure)
    return localization_matrix(sp.monomials, sp.polynomial, moment_vars)
end

function primal_scheme_constraint(model::JuMP.Model, sp::SchemePart, moment_vars::MM.Measure)
    cref = @constraint model primal_scheme_part(sp, moment_vars) in sp.moi_set
    return cref
end

export PutinarScheme
"""
PutinarScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)

Return a Putinar Approximation Scheme.
"""
struct PutinarScheme <: AbstractApproximationScheme
    equality::AbstractEqualityHandler
    sparsity::SumOfSquares.Sparsity
    basis_type::Type{<:MB.AbstractPolynomialBasis}
    function PutinarScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)
        new(equality, sparsity, basis)
    end
end

# this is very inefficient but does what I want.
function approximation_scheme(scheme::PutinarScheme, K::AbstractBasicSemialgebraicSet, vars::Vector{<: MP.AbstractVariable}, d::Int)

    ineqs = [prod(vars.^0)]

    if !isempty(inequalities(K))
        ineqs = [ineqs..., inequalities(K)...]
    end

    if scheme.equality isa SplitEquality
        for eq in equalities(K)
            append!(ineqs, [eq, -eq])
            eqs = typeof(ineqs)()
        end
    elseif scheme.equality isa ReduceEquality
        @error "ask Benoit"
    else
        eqs = equalities(K)
    end

    schemeparts = SchemePart[]

    for ineq in ineqs
        deg = Int(floor((d-maxdegree(ineq))/2))
        mons = maxdegree_basis(scheme.basis_type, vars, deg)
        if length(mons) == 1 
            moiset = MOI.Nonnegatives(1)
        else
            moiset = PSDCone()
        end
        push!(schemeparts, SchemePart(ineq, mons, moiset))
    end
    for eq in eqs
        deg = d-maxdegree(eq)
        mons = maxdegree_basis(scheme.basis_type, vars, deg)
        push!(schemeparts, SchemePart(eq, mons, MOI.Zeros(length(mons))))
    end
    return schemeparts
end

export SchmuedgenScheme
"""
SchmuedgenScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)

Return a Schmuedgen Approximation Scheme.

"""
struct SchmuedgenScheme <: AbstractApproximationScheme
    equality::AbstractEqualityHandler
    sparsity::SumOfSquares.Sparsity
    basis_type::Type{<:MB.AbstractPolynomialBasis}
    function SchmuedgenScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)
        new(equality, sparsity, basis)
    end
end

function approximation_scheme(scheme::SchmuedgenScheme, K::AbstractBasicSemialgebraicSet, vars::Vector{<: MP.AbstractVariable}, d::Int)

    ineqs = [prod(vars.^0)]

    if !isempty(inequalities(K))
        ineqs = [ineqs..., inequalities(K)...]
    end

    if scheme.equality isa SplitEquality
        for eq in equalities(K)
            append!(ineqs, [eq, -eq])
            eqs = typeof(ineqs)()
        end
    elseif scheme.equality isa ReduceEquality
        @error "ask Benoit"
    else
        eqs = equalities(K)
    end

    schemeparts = SchemePart[]
    i = 1
        ineq = ineqs[i]
        deg = Int(floor((d-maxdegree(ineq))/2))
        mons = maxdegree_basis(scheme.basis_type, vars, deg) 
        if length(mons) == 1
            moiset = MOI.Nonnegatives(1)
        else
            moiset = PSDCone()
        end
        push!(schemeparts, SchemePart(ineq, mons, moiset))

    for i = 2:length(ineqs)
        ineq = ineqs[i]
        deg = Int(floor((d-maxdegree(ineq))/2))
        mons = maxdegree_basis(scheme.basis_type, vars, deg) 
        if length(mons) == 1
            moiset = MOI.Nonnegatives(1)
        else
            moiset = PSDCone()
        end

        push!(schemeparts, SchemePart(ineq, mons, moiset))
        for j = i+1:length(ineqs)
            if maxdegree(ineqs[i]) + maxdegree(ineqs[j]) <= d
                ineq = ineqs[i]*ineqs[j]
                deg = Int(floor((d-maxdegree(ineq))/2))
                mons = maxdegree_basis(scheme.basis_type, vars, deg)
                if length(mons) == 1
                    moiset = MOI.Nonnegatives(1)
                else
                    moiset = PSDCone()
                end
                push!(schemeparts, SchemePart(ineq, mons, moiset))
            end
        end
    end

    for eq in eqs
        deg = d-maxdegree(eq)
        mons = maxdegree_basis(scheme.basis_type, vars, deg)
        push!(schemeparts, SchemePart(eq, mons, MOI.Zeros(length(mons))))
    end
    return schemeparts
end

