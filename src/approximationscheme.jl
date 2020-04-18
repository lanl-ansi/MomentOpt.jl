export KeepEquality, ReduceEquality, SplitEquality
abstract type AbstractEqualityHandler end
struct KeepEquality <: AbstractEqualityHandler end
struct ReduceEquality <: AbstractEqualityHandler end
struct SplitEquality <: AbstractEqualityHandler end


export PutinarScheme
"""
    PutinarScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)

Return a Putinar Approximation Scheme.
"""
struct PutinarScheme <: AbstractApproximationScheme
    equality::AbstractEqualityHandler
    sparsity::SumOfSquares.Sparsity
    basis_type #todo figure type/ask Benoit
    function PutinarScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)
        new(equality, sparsity, basis)
    end
end

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
    
    dict = Dict{Any, Any}()
    for ineq in ineqs
        deg = Int(floor((d-maxdegree(ineq))/2))
        mons = monomials(maxdegree_basis(scheme.basis_type, vars, deg))
        if length(mons) == 1
            moiset = MOI.GreaterThan(0.0)
            mons = mons
        else
            moiset = PSDCone()
        end
        dict[ineq*mons*transpose(mons)] = moiset
    end
    for eq in eqs
        deg = d-maxdegree(eq)
        dict[monomials(maxdegree_basis(scheme.basis_type, vars, deg)).*eq] = MOI.GreaterThan(0)
    end
    return dict
end

export SchmuedgenScheme
"""
    SchmuedgenScheme(;equality = KeepEquality(), sparsity = NoSparsity(), basis = MonomialBasis)

Return a Schmuedgen Approximation Scheme.

"""
struct SchmuedgenScheme <: AbstractApproximationScheme
    equality::AbstractEqualityHandler
    sparsity::SumOfSquares.Sparsity
    basis_type #todo figure type/ask Benoit
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
    dict = Dict{Any, Any}()
    for i = 1:length(ineqs)
        ineq = ineqs[i]
        deg = Int(floor((d-maxdegree(ineq))/2))
        mons = monomials(maxdegree_basis(scheme.basis_type, vars, deg)) 
        if length(mons) ==1
            moiset = MOI.GreaterThan(0.0)
            mons = first(mons)
        else
            moiset = PSDCone()
        end
        dict[[ineq.*mons*transpose(mons)]] = moiset
        for j = i+1:length(ineqs)
            if maxdegree(ineqs[i]) + maxdegree(ineqs[j]) <= d
                ineq = ineqs[i]*ineqs[j]
                deg = Int(floor((d-maxdegree(ineq))/2))
                mons = monomials(maxdegree_basis(scheme.basis_type, vars, deg))
                if length(mons) ==1
                    moiset = MOI.GreaterThan(0.0)
                    mons = first(mons)
                else
                    moiset = PSDCone()
                end
                dict[ineq.*mons*transpose(mons)] = moiset
            end
        end
    end
    for eq in eqs
        deg = d-maxdegree(eq)
        dict[monomials(maxdegree_basis(scheme.basis_type, vars, deg)).*eq] = MOI.GreaterThan(0)
    end
    return dict
end

