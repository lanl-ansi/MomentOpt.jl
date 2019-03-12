#created: TW 2019/02/22
#modified: TW 2019/02/27
#
#

#using Mosek
using SDPA

function shcl_01(pde,relaxation_order)

t_min = 0 # may not be changed
t_max = pde[:bounds][:t_max]
x_min = pde[:bounds][:x_min]
x_max = pde[:bounds][:x_max]
y_min = pde[:bounds][:y_min]
y_max = pde[:bounds][:y_max]

m = MomentModel(solver=SDPASolver())

@polyvar t
@polyvar x
@polyvar y
@polyvar v

scaling = t_max/(x_max-x_min)
println(scaling)

# functions for Lebesque marginals
lebt_mom(i) = 1/(i+1)
lebx_mom(j) = ((1/2)^(j+1)-(-1/2)^(j+1))/(j+1)
lebtx_mom(i,j) = lebt_mom(i)*lebx_mom(j) 


# occupation measure
add_measure!(m,:mu)
add_support!(m, :mu, @set((1.0-t)*t>=0 && (1.0/2-x)*(x+1.0/2)>=0 && (y_max-y)*(y-y_min)>=0)) 
mom_fun = (i,j)->t^i*x^j
add_mfun_constraint!(m,[(mom_fun,:mu)],:eq,lebtx_mom,2) 

# initial measure 
add_measure!(m,:mu_init)
add_support!(m,:mu_init, @set(t==0 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0))
mom_funy = (i,j,k)-> t^i*x^j*y^k
add_mfun_constraint!(m,[(mom_funy,:mu_init)],:eq,(i,j,k)->0^i*quadgk(xx-> xx^j*pde[:initial](xx)^k,-1/2,1/2)[1],3) #TODO: remove measure from formulation)


# terminal measure 
add_measure!(m,:mu_term)
add_support!(m, :mu_term, @set(t==1  && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0)) 
add_mfun_constraint!(m, [(mom_fun,:mu_term)], :eq, (i,j) -> 1^i*lebx_mom(j) ,2)


# left boundary
add_measure!(m,:mu_left)
add_support!(m,:mu_left, @set(t*(1-t)>=0 && x==-1/2 && (y_max-y)*(y-y_min)>=0))
if haskey(pde[:boundaries],:left)
	add_mfun_constraint!(m,[(mom_funy,:mu_left)],:eq,(i,j,k)->(-1/2)^j*quadgk(tt->tt^i*pde[:boundaries][:left](tt)^k,0,1)[1],3) 
else
	add_mfun_constraint!(m,[(mom_fun,:mu_left)],:eq,(i,j)->lebt_mom(i)*(-1/2)^j,2)
end

# right boundary
add_measure!(m,:mu_right)
add_support!(m,:mu_right, @set(t*(1-t)>=0 && x==1/2 && (y_max-y)*(y-y_min)>=0))
if haskey(pde[:boundaries],:right)
	add_mfun_constraint!(m,[(mom_funy,:mu_right)],:eq,(i,j,k)->(1/2)^j*quadgk(tt->tt^i*pde[:boundaries][:right](tt)^k,0,1)[1],3)
else
	add_mfun_constraint!(m,[(mom_fun,:mu_right)],:eq,(i,j)->lebt_mom(i)*(1/2)^j,2)
end

# periodic condition
if haskey(pde[:boundaries],:periodic) && pde[:boundaries][:periodic]
	mom_fun_nox = (i,k) -> t^i*y^k
	neg_mom_fun_nox = (i,k) -> -t^i*y^k
	add_mfun_constraint!(m,[(mom_fun_nox,:mu_right),(neg_mom_fun_nox,:mu_left)],:eq,(i,k)->0,2)
end


# conservation law
phi_con_mu = (i,j) ->  -(differentiate(t^i*x^j,t)*y+differentiate(t^i*x^j,x)*scaling*pde[:flux](y))
momy = (i,j) -> t^i*x^j*y
mmomy = (i,j) -> -momy(i,j)
momfy = (i,j) -> scaling*t^i*x^j*pde[:flux](y)
mmomfy = (i,j) -> -momfy(i,j)
add_mfun_constraint!(m,[(momy,:mu_term),(mmomy,:mu_init),(momfy,:mu_right),(mmomfy,:mu_left),(phi_con_mu,:mu)],:eq,(i,j)->0,2)

# Preparation to use Kruzkhov

mom_funv = (i,j,k,l) -> t^i*x^j*y^k*v^l
mom_funV = (i,j,k,l) -> -t^i*x^j*y^k*(y_max^(l+1)-y_min^(l+1))/(l+1)/(y_max-y_min)

# auxiliary measures for sign(y-v) 
add_measure!(m,:nu_p)
add_support!(m,:nu_p, @set(t*(1-t)>=0 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && y-v>=0 )) 
add_measure!(m,:nu_m)
add_support!(m,:nu_m, @set(t*(1-t)>=0 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && v-y>=0 ))
add_mfun_constraint!(m, [(mom_funv,:nu_p),(mom_funv,:nu_m),(mom_funV,:mu)],:eq,(i,j,k,l)->0,4)

add_measure!(m,:nu_init_p)
add_support!(m,:nu_init_p, @set(t==0 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && y-v>=0 )) 
add_measure!(m,:nu_init_m)
add_support!(m,:nu_init_m, @set(t==0 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && v-y>=0 ))
add_mfun_constraint!(m, [(mom_funv,:nu_init_p),(mom_funv,:nu_init_m),(mom_funV,:mu_init)],:eq,(i,j,k,l)->0,4)

add_measure!(m,:nu_term_p)
add_support!(m,:nu_term_p, @set(t==1 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && y-v>=0 )) 
add_measure!(m,:nu_term_m)
add_support!(m,:nu_term_m, @set(t==1 && (1/2-x)*(x+1/2)>=0 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && v-y>=0 ))
add_mfun_constraint!(m, [(mom_funv,:nu_term_p),(mom_funv,:nu_term_m),(mom_funV,:mu_term)],:eq,(i,j,k,l)->0,4)

add_measure!(m,:nu_right_p)
add_support!(m,:nu_right_p, @set(t*(1-t)>=0 && x==1/2 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && y-v>=0 )) 
add_measure!(m,:nu_right_m)
add_support!(m,:nu_right_m, @set(t*(1-t)>=0 && x==1/2 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && v-y>=0 ))
add_mfun_constraint!(m, [(mom_funv,:nu_right_p),(mom_funv,:nu_right_m),(mom_funV,:mu_right)],:eq,(i,j,k,l)->0,4)

add_measure!(m,:nu_left_p)
add_support!(m,:nu_left_p, @set(t*(1-t)>=0 && x==-1/2 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && y-v>=0 )) 
add_measure!(m,:nu_left_m)
add_support!(m,:nu_left_m, @set(t*(1-t)>=0 && x==-1/2 && (y_max-y)*(y-y_min)>=0 && (y_max-v)*(v-y_min)>=0 && v-y>=0 ))
add_mfun_constraint!(m, [(mom_funv,:nu_left_p),(mom_funv,:nu_left_m),(mom_funV,:mu_left)],:eq,(i,j,k,l)->0,4)


# Kruzkhov entropies
add_measure!(m,:nu_kruzk) 
add_support!(m,:nu_kruzk, @set( t*(1-t)>=0 && (1/2-x)*(x+1/2)>=0 && (y_max-v)*(v-y_min)>=0 ))

eta_p = (i,j,l) ->  t^i*x^j*v^l*(y-v)
eta_m = (i,j,l) ->  t^i*x^j*v^l*(v-y)
q_p = (i,j,l) ->  t^i*x^j*v^l*scaling*(pde[:flux](y)-pde[:flux](v))
q_m = (i,j,l) ->  t^i*x^j*v^l*scaling*(pde[:flux](v)-pde[:flux](y))
eta_p_neg = (i,j,l) -> -eta_p(i,j,l)
eta_m_neg = (i,j,l) -> -eta_m(i,j,l)
q_p_neg = (i,j,l) ->  -q_p(i,j,l)
q_m_neg = (i,j,l) ->  -q_m(i,j,l)
d_p = (i,j,l) -> -(differentiate(eta_p(i,j,l),t)+differentiate(q_p(i,j,l),x))
d_m = (i,j,l) -> -(differentiate(eta_m(i,j,l),t)+differentiate(q_m(i,j,l),x))
mom_fun_kruzk = (i,j,l) -> t^i*x^j*v^l
add_mfun_constraint!(m,[(d_p,:nu_p),(d_m,:nu_m),(eta_p_neg,:nu_init_p),(eta_m_neg,:nu_init_m),(eta_p,:nu_term_p),(eta_m,:nu_term_m),(q_p_neg,:nu_left_p),(q_m_neg,:nu_left_m),(q_p,:nu_right_p),(q_m,:nu_right_m),(mom_fun_kruzk,:nu_kruzk)],:eq,(i,j,l)->0,3)


# objective
foo = sum([t^(2*i) + x^(2*i)+ y^(2*i) for i=1:relaxation_order])
add_objective!(m,:Min,[(foo,:mu)])

mm_solve(m,relaxation_order)

return m
end

