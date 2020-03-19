export @measure

function _macro_assign_and_return(code, variable, name, model_for_registering)
    #TODO look at JuMP version and extend.
    return quote
        $variable = $code
        $(esc(name)) = $variable
    end
end


macro measure(args...)
    # adapted from JuMP.@variable
   
    # Initiate an error function
    _error(str...) = JuMP._macro_error(:measure, args, str...)

    # The first argument is always assumed to be the GMPModel the measure should be added to.
    # Therefore we need to escape the first element of the inputs (Vector of Expr).
    # The result model now is an Expression that - evaluate during the macro evaluation -
    # refers to the model in the scope where the macro is called, i.e., to the thing the user 
    # has entered first after @measure
    model = esc(args[1]) 

    # sorting the remaining arguments into extra, kw_args and requested container
    extra, kw_args, requestedcontainer = JuMP.Containers._extract_kw_args(args[2:end])

    # we need to generate at least two inputs in extra
    # the first element is interpretate as the name
    # the second one is interpretate as the variables
    # if there is only a single element in extra this is assumed to be 
    # the polynomial variables and the measure variable is anonymous
    
    if isempty(extra)
        return :(throw(ArgumentError("Specify at least model and polynomial variables.\nE.g., @measure model variables")))
    end
              
    if length(extra) == 1
        var = gensym()
        anon_singleton = true
        polyvars = extra[1]
    else
        var = popfirst!(extra)
        anon_singleton = false
        polyvars = popfirst!(extra)
    end
    
    # Similar to the model, the polyvars expression should refer to the vector of polynomial
    # variables which was entered in the macro call. We therfore need to escape it too.
    polyvars = esc(polyvars)

   # @info var
   # @info anon_singleton
   # @info polyvars

    # Sorting out measure key words
    info_kw_args = filter(_is_measure_info_keyword, kw_args)

    # Sort out base_name keyword for containers
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    @assert length(base_name_kw_args)<=1 "Cannot specify two different base_names."

    # Building the info expression for new variable
    # JuMP._keywordify escapes the second argument of the keywords
    infoexpr = _MeasureInfoExpr(polyvars; JuMP._keywordify.(info_kw_args)...)

   # @info infoexpr
    
    # There are four cases to consider:
    # x                                       | type of x | x.head
    # ----------------------------------------+-----------+------------
    # var                                     | Symbol    | NA
    # var[1:2]                                | Expr      | :ref


    anonvar = JuMP.isexpr(var, :vect) || JuMP.isexpr(var, :vcat) || anon_singleton
    variable = gensym()
    name = JuMP.Containers._get_name(var)
    # @info name
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        Base.error("Expression $name should not be used as a variable name. Use the \"anonymous\" syntax $name = @measure(model, ...) instead.")
    end

    info = JuMP._constructor_expr(infoexpr)

    # @info info
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
    buildcall = :( build_variable($_error, $info))    
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
        macro_code = _macro_assign_and_return(creation_code, variable, name, model)
    end

    return macro_code
end


