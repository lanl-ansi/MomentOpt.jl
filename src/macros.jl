export @measure, @mconstraint, @mconstraints, @mobjective


macro measure(gmp,mu)
	Expr(:block,:(add_measure!($(esc(gmp)),$(esc(mu)))))
end


macro measure(gmp, name, vars)
	quote
        	$(esc(name))= Measure($(string(name)),$(esc(vars)))
		add_measure!($(esc(gmp)),$(esc(name)))
	end
end

macro measure(gmp, name, vars, kwargs...)
	quote
        	$(esc(name))= Measure($(string(name)),$(esc(vars));$(esc(kwargs...)))
		add_measure!($(esc(gmp)),$(esc(name)))
	end
end

macro mconstraint(gmp,momcons)
	quote
		add_constraint!($(esc(gmp)),$(esc(momcons)))
	end
end

macro mconstraints(gmp,momcons)
	quote
		add_constraints!($(esc(gmp)),$(esc(momcons)))
	end
end

macro mobjective(gmp,args...)
	cargs = collect(args)
	nargs = length(cargs)
	if nargs==1
	return quote
		add_objective!($(esc(gmp)),$(esc(cargs[1])))
	end
	elseif nargs==2
	return quote
		add_objective!($(esc(gmp)),$(esc(cargs[1])),$(esc(cargs[2])))
	end
	else
	return quote
		add_objective!($(esc(gmp)),$(esc(cargs[1])),$(esc(cargs[2])),$(esc(cargs[3])))
	end
	end
end
