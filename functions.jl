using MultivariatePolynomials
using JuMP
using PolyJuMP
using SumOfSquares
using DynamicPolynomials


function MomentModel(;solver = solvertag)
    mm = Dict()
    mm[:meas_list] = []
    mm[:supp_list] = Dict()
    mm[:objective] = Dict()
    mm[:cons_list] = []
    mm[:solver] = solver
    return mm
end

function add_measure!(mmodel, meas)
    	push!(mmodel[:meas_list],meas)
    return mmodel
end

function rm_measure!(mmodel, meas)
    pop!(mmodel[:meas_list],meas)
end

function add_support!(mmodel, meas, saset)
    mmodel[:supp_list][meas]= saset
    return mmodel
end

function add_objective!(mmodel, mode, tuple_list)
    # mode is :Min, or :Max
    mmodel[:objective][:mode] = mode
    for (pol,meas) in tuple_list
        mmodel[:objective][meas] = pol
    end
    return mmodel
end

function add_mconstraint!(mmodel, dict, mode, RHS)
    # mode is :leq, :geq, or :eq
    constr = Dict()
    constr[:type] = :fin
    constr[:cons] = Dict()
    for (polvec,meas) in dict
        constr[:cons][meas]=polvec
    end
    constr[:mode] = mode
    constr[:RHS] = RHS

    push!(mmodel[:cons_list],constr)
    return mmodel
end

function add_mfun_constraint!(mmodel, int_tuple2pol_list, mode, int_tuple2real, int_tuple_size)
    # mode is :leq, :geq, or :eq
    constr = Dict()
    constr[:type] = :inf
    constr[:mode] = mode
    constr[:tuple_size] = int_tuple_size
    constr[:RHS] = int_tuple2real
    constr[:cons] = Dict()
    for (fun,meas) in int_tuple2pol_list
        constr[:cons][meas] = fun
    end
    push!(mmodel[:cons_list],constr)
    return mmodel
end

function hompow(n,d)
    if n > 1
	if d>0
	    ct=0
	    AUX = zeros(Int,binomial(n+d-1,d),n)
	    for d0 = d:-1:0
		aux = hompow(n-1,d-d0)
		aux_s = size(aux,1)
		AUX[ct+1:ct+aux_s,1]=Int(d0)*ones(Int,aux_s,1)
		AUX[ct+1:ct+aux_s,2:n]=aux
		ct = ct + aux_s
	    end
	    return AUX
	else
	    return zeros(Int,1,n)
	end
    else
	return [Int(d)]
    end 
end


function truncate_constraint(mmodel,constr, order, maxct)
    # very bad function!
    # mode is :leq, :geq, or :eq
    tconstr = Dict()
    tconstr[:type] = :fin
    tconstr[:mode] = constr[:mode]
    tconstr[:RHS] = []
    tconstr[:cons] = Dict()

    for meas in mmodel[:meas_list]
        if haskey(constr[:cons],meas)
            tconstr[:cons][meas] = []
        end
    end
    cct = 0;
    while cct<=maxct
        bool = true
	int_tuples = hompow(constr[:tuple_size],cct)
	s_int_tuples = size(int_tuples,1)
	k=0
	while(k<s_int_tuples)&&bool
	    k= k+1
            for meas in mmodel[:meas_list]
            	if haskey(constr[:cons],meas)
		    bool = bool && (maxdegree(constr[:cons][meas](int_tuples[k,:]...))<=2*order)
            	end
            end
            if bool
                for meas in mmodel[:meas_list]
                    if haskey(constr[:cons],meas)
                	append!(tconstr[:cons][meas],[constr[:cons][meas](int_tuples[k,:]...)])
               	    end
            	end
            	append!(tconstr[:RHS],constr[:RHS](int_tuples[k,:]...))
            else
            	cct = maxct+1
            end
	end
        cct=cct+1
    end
    return tconstr
end

#=
function truncate_constraint(mmodel,constr, order, maxct)
    # very bad function!
    # mode is :leq, :geq, or :eq
    tconstr = Dict()
    tconstr[:type] = :fin
    tconstr[:mode] = constr[:mode]
    tconstr[:RHS] = []
    tconstr[:cons] = Dict()

    cct = 1;
    for meas in mmodel[:meas_list]
        if haskey(constr[:cons],meas)
            tconstr[:cons][meas] = []
        end
    end

    while cct<=maxct
        bool = true
        for meas in mmodel[:meas_list]
            if haskey(constr[:cons],meas)
                bool = bool && (maxdegree(constr[:cons][meas](cct))<=2*order)
            end
        end
        if bool
            for meas in mmodel[:meas_list]
                if haskey(constr[:cons],meas)
                    append!(tconstr[:cons][meas],[constr[:cons][meas](cct)])
                end
            end
            append!(tconstr[:RHS],constr[:RHS](cct))
            cct=cct+1
        else
            cct = maxct+1
        end
    end
    return tconstr
end
=#

function check_constraint(mmodel,constr, order)
    bool = true
    for meas in mmodel[:meas_list]
	    if haskey(constr[:cons],meas)
	        for i =1:length(constr[:RHS])
                bool = bool && (maxdegree(constr[:cons][meas][i])<=2*order)
		end
        end
    end
    if bool
        return constr
    else
        println("Increase relaxation order!")
        return
    end
end


function get_constraints(mmodel, order)
    constraints = Dict()
    constraints[:RHS] = []
    constraints[:mode] = []
    for meas in mmodel[:meas_list]
        constraints[meas]=[]
    end
    for constr in mmodel[:cons_list]
        if constr[:type] == :inf
            tconstr = truncate_constraint(mmodel,constr, order,1000)
        else
            tconstr = check_constraint(mmodel,constr,order)
        end
        append!(constraints[:RHS], tconstr[:RHS])
        append!(constraints[:mode], [tconstr[:mode] for i=1:length(tconstr[:RHS])])
        n_cons = length(tconstr[:RHS])
        for meas in mmodel[:meas_list]
            # append!(constraints[meas],get(tconstr[:cons],meas,zeros(Polynomial{true,Float64}, n_cons)))
	      append!(constraints[meas],get(tconstr[:cons],meas,zeros(n_cons)))
        end
    end
# should give a warning when degree of constraints is higher than requested order
# println(constraints)
    return constraints
end

function mm_solve(mmodel,order)
    bool = true
    for meas in mmodel[:meas_list]
        for sup_pol in mmodel[:supp_list][meas].p
            bool = bool && maxdegree(sup_pol)<=2*order
        end
        bool = bool && maxdegree(get(mmodel[:objective],meas,convert(Polynomial{true,Float64},0)))<=2*order
    end

    if ~bool
        println("Increase relaxation order!")
        return
    end
    constraints = get_constraints(mmodel, order)
    n_cons = length(constraints[:RHS])
    m = SOSModel(solver = mmodel[:solver])

    p = Dict()
    for i = 1:n_cons
        if  constraints[:mode][i] == :eq
            p[i] = @variable(m)
        elseif ((constraints[:mode][i] == :leq)&&(mmodel[:objective][:mode]==:Max))||((constraints[:mode][i] == :geq)&&(mmodel[:objective][:mode]==:Min))
            p[i] = @variable(m, lowerbound = 0)
        else
            p[i] = @variable(m, upperbound = 0)
        end
        # doublecheck ten times if the inequalities for p are right
    end
    mmodel[:dual_list] = Dict()
    if  (mmodel[:objective][:mode]==:Max)
        @objective m Min sum(p[i]*constraints[:RHS][i] for i in 1:n_cons)
        for meas in mmodel[:meas_list]
            expr = sum(p[i]*constraints[meas][i] for i in 1:n_cons)-get(mmodel[:objective],meas,convert(Polynomial{true,Float64},0))
            supp = mmodel[:supp_list][meas]
            mmodel[:dual_list][meas] = PolyJuMP.addpolyconstraint!(m, expr, SOSCone(), supp , MonomialBasis; maxdegree = 2*order)
        end
    else
        @objective m Max sum(p[i]*constraints[:RHS][i] for i in 1:n_cons)
        for meas in mmodel[:meas_list]
            expr = get(mmodel[:objective],meas,convert(Polynomial{true,Float64},0))-sum(p[i]*constraints[meas][i] for i in 1:n_cons)
            supp = mmodel[:supp_list][meas]
            mmodel[:dual_list][meas] = PolyJuMP.addpolyconstraint!(m, expr , SOSCone(), supp, MonomialBasis; maxdegree = 2*order)
        end
    end

    status = solve(m)
    mmodel[:sosmodel] = m
    mmodel[:moments] = Dict()
    for meas in mmodel[:meas_list]
        mmodel[:moments][meas] = get_moments(mmodel,meas)
    end
    return status
end

function get_moments(mmodel, meas)
    sos_moments = getdual(mmodel[:dual_list][meas])
    moments = Dict()
    moments[:basis] = sos_moments.x
    if mmodel[:objective][:mode] == :Min
        moments[:moments] = -sos_moments.a
    else
        moments[:moments] =  sos_moments.a
    end
    return moments
end
