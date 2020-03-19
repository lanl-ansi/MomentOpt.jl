export @measure

function JuMP._assert_valid_model(m, macrocode)
    # assumes m is already escaped
    quote
        JuMP._valid_model($m, $(JuMP.quot(m.args[1])))
        $macrocode
    end
end

function JuMP._macro_assign_and_return(code, variable, name;
                                  model_for_registering=nothing)
    return quote
        $(if model_for_registering !== nothing
            :(_error_if_cannot_register($model_for_registering,
                                        $(quot(name))))
          end)
        $variable = $code
        $(if model_for_registering !== nothing
            :(object_dictionary($model_for_registering)[$(quot(name))] =
              $variable)
          end)
        # This assignment should be in the scope calling the macro
        $(esc(name)) = $variable
    end
end


macro measure(args...)
    _error(str...) = JuMP._macro_error(:measure, args, str...)

    model = esc(args[1])

    extra, kw_args, requestedcontainer = JuMP.Containers._extract_kw_args(args[2:end])

    # we need to generate at least two inputs in extra
    # the first element is interpretate as the name
    # the second one is interpretate as the variables
    # if there is only a single element in extra this is assumed to be 
    # the polynomial variables and the measure variable is anonymous
    
    @assert !isempty(extra)

    if length(extra) == 1
        var = gensym()
        anon_singleton = true
        polyvars = extra[1]
    else
        var = popfirst!(extra)
        anon_singleton = false
        polyvars = popfirst!(extra)
    end

    info_kw_args = filter(_is_measure_info_keyword, kw_args)
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    infoexpr = _MeasureInfoExpr(; JuMP._keywordify.(info_kw_args)...)

    # There are four cases to consider:
    # x                                       | type of x | x.head
    # ----------------------------------------+-----------+------------
    # var                                     | Symbol    | NA
    # var[1:2]                                | Expr      | :ref


    anonvar = JuMP.isexpr(var, :vect) || JuMP.isexpr(var, :vcat) || anon_singleton
    variable = gensym()
    name = JuMP.Containers._get_name(var)
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        Base.error("Expression $name should not be used as a variable name. Use the \"anonymous\" syntax $name = @measure(model, ...) instead.")
    end

    info = JuMP._constructor_expr(infoexpr)
    if isa(var, Symbol)
        # Easy case - a single variable
        name_code = base_name
    else
        isa(var, Expr) || _error("Expected $var to be a variable name")
        # We now build the code to generate the variables (and possibly the
        # SparseAxisArray to contain them)
        idxvars, indices = JuMP.Containers._build_ref_sets(_error, var)
        name_code = JuMP._name_call(base_name, idxvars)
    end

    # Code to be used to create each variable of the container.
    buildcall = :( build_variable($_error, $info, $(extra...)) )
    variablecall = :( add_variable($model, $buildcall, $name_code) )
    
    if isa(var, Symbol)         
        # The looped code is trivial here since there is a single variable
        creation_code = variablecall
    else
        creation_code = JuMP.Containers.container_code(idxvars, indices, variablecall, requestedcontainer)
    end

    if anonvar
        # Anonymous variable, no need to register it in the model-level
        # dictionary nor to assign it to a variable in the user scope.
        # We simply return the variable
        macro_code = creation_code
    else
        # We register the variable reference to its name and
        # we assign it to a variable in the local scope of this name
        macro_code = JuMP._macro_assign_and_return(creation_code, variable, name,
                                              model_for_registering = model)
    end

    println(creation_code)
    return JuMP._assert_valid_model(model, macro_code)
end
