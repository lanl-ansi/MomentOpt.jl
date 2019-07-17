export @measure, @mconstraint, @mconstraints


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

