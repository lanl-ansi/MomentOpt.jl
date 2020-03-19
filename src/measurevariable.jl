export MeasureVariable, Measure

abstract type AbstractRelaxationType end
struct Putinar <: AbstractRelaxationType end #TODO move

mutable struct _MeasureInfoExpr
    poly_variables::Any
    support::Any
    moment_basis::Any
    relax_type::Any
end

function _is_measure_info_keyword(kw::Expr)
    kw.args[1] in [:support, :momentbasis, :relaxtype]
end

function _MeasureInfoExpr(;
                          polyvars = [], 
                          support = FullSpace(),
                          momentbasis = MB.MonomialBasis,
                          relaxtype = Putinar)
    return _MeasureInfoExpr(polyvars, support, momentbasis, relaxtype)
end

mutable struct MeasureInfo{PT, SAS, BT, ART}
    poly_variables::Vector{PT}
    support::SAS
    moment_basis::BT
    relax_type::ART
end

function JuMP._constructor_expr(info::_MeasureInfoExpr)
    return :(MeasureInfo($(info.poly_variables), 
                         $(info.support), 
                         $(info.moment_basis),
                         $(info.relax_type)
                        )
            )
end

mutable struct MeasureVariable <: JuMP.AbstractVariable
    info::MeasureInfo
end

function JuMP.build_variable(_error::Function, info::MeasureInfo; extra_kw_args...)
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    return MeasureVariable(info)
end

const Measure = MeasureVariable
