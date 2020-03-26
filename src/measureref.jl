export GMPRef

struct GMPRef <: JuMP.AbstractVariableRef
    model::GMPModel
    index::Int
end

function GMPRef(m::GMPModel)
    index = MOI.add_variable(model_info(m))
    return GMPRef(m, index)
end

function JuMP.add_variable(model::GMPModel, v <: AbstractGMPType, name::String="")
    var_ref = GMPRef(model)
    model_data(model).variables[var_ref] = v
    if !isempty(name)
        set_name(var_ref, name)
    end
    return var_ref   
end
function JuMP.variable_by_name(model::GMPModel, name::String)
    index = findfirst(x -> x == name, model_info(model).variable_names)
    if index isa Nothing
        return nothing
    else
        return GMPRef(model, index)
    end
end
