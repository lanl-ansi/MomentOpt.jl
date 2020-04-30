"""
 Approximate the volume of a semi algebraic set with MomentOpt.jl.
 For theoretical background see : https://homepages.laas.fr/henrion/Papers/volume.pdf

 Let K ⊆ B ⊆ R^n be semialgebraic compact sets and the moments of a the Lebesgue meausre 
 on B be known. One is interested in computing the volume of K. This can be done by solving
 the Generalized Moment Problem

  max <μ,1>
 	μ + ν = λ
	μ ∈ M(K)
	ν ∈ M(B)
 where λ is the Lebesgue measure on B.

 The equality of measures can be understood in the weak sense, i.e,
	<μ,φ> + <ν,φ> = <λ,φ>
 for all test functions φ. 
 As B is compact one can choose the polynomials, or any basis of the space of polynomials as
 test functions. In the code below we choose the monomials as test functions. 

 When choosing the monomials up to a fixed degree as test functions, the dual problem reads

 min ∑p_α <λ,x^α>
     ∑p_αx^α ⩾ 1 on K
     ∑p_αx^α ⩾ 0 on B 

 i.e., the optimal polynomial p = ∑p_αx^α is an over approximation of the indicator function of 
 K on B.
"""

using DynamicPolynomials
using MomentOpt
using MosekTools


# relaxation order for the Generalized Moment Problem
order = 10

# polnomial variables
@polyvar x y

# define set K of which the volume should be computed
K = @set(1-x^2-y^2>=0)
# define superset of K where the moments of the Lebegue measure can be computed
B = @set(1-x^2>=0 && 1-y^2>=0)

# define measures on K and B
gmp = GMPModel()

@variable gmp μ Meas([x,y], support = K)
@variable gmp ν Meas([x,y], support = B)

# We want to maximize the mass of mu
@objective gmp Max Mom(1, μ)

# The constraint is mu + nu = Lebesgue on B
λ = lebesgue_measure_box(variable_box(x => [-1, 1], y => [-1, 1]); normalize = true)
con = @constraint gmp μ + ν == λ

# In order to relax the problem, we need to specify the relaxation order and a solver factory:
set_approximation_mode(gmp, DUAL_STRENGTHEN_MODE())
set_optimizer(gmp, Mosek.Optimizer)
set_approximation_degree(gmp, 2*order)
optimize!(gmp)

# The optimal value is an approximation to the volume of K (which is π in this case)
println("Volume of approximation: $(objective_value(gmp)*4)")
println(termination_status(gmp))
#=
# The polynomial over approximation of the indicator function can be defined
# from the dual solution as follows:
poly = MomentOpt.approx_crefs(gmp)[index(con)].dual
println()
println(
"""
The approximated function can be plotted, e.g., using the following commands:

using Plots
pyplot()
xx = yy = range(-1, stop = 1, length = 100)
f(xx,yy) = poly(x=>xx,y=>yy)
plot(xx, yy, f, st= :surface)
""")

# The convergence can be improved by adding "Stokes" constraints.
mons = monomials([x,y],0:2*order-1)
#pons = polynomial.(mons)

@constraint gmp stokes_x[i=1:length(mons)] Mom(differentiate(mons[i]*(1-x^2-y^2), x), μ) == 0
@constraint gmp stokes_y[i=1:length(mons)] Mom(differentiate(mons[i]*(1-x^2-y^2), y), μ) == 0

optimize!(gmp)
println("Approximation with Stokes: $(objective_value(gmp)*4)")
println(termination_status(gmp))
=#
