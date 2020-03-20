export MeasureVariable, Measure 

abstract type AbstractRelaxationType end
struct Putinar <: AbstractRelaxationType end #TODO move and change
struct Handelman <: AbstractRelaxationType end #TODO move and change

mutable struct _MeasureInfoExpr
    poly_variables::Any
    support::Any
    moment_basis::Any
    relax_type::Any
end

function _is_measure_info_keyword(kw::Expr)
    kw.args[1] in [:support, :moment_basis, :relax_type]
end

function _MeasureInfoExpr(polyvars;
                          support = :(FullSpace()),
                          moment_basis = :(MB.MonomialBasis),
                          relax_type = :(Putinar()))
        _MeasureInfoExpr(polyvars, support, moment_basis, relax_type)
end

mutable struct MeasureInfo
    poly_variables
    support
    moment_basis
    relax_type
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

info(mv::MeasureVariable) = mv.info
poly_variables(mv::MeasureVariable) = info(mv).poly_variables

function JuMP.build_variable(_error::Function, info::MeasureInfo; extra_kw_args...)
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    return MeasureVariable(info)
end

const Measure = MeasureVariable
