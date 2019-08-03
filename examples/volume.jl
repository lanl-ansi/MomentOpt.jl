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
using SemialgebraicSets

using MomentOpt

using MosekTools
using Plots
pyplot()

# relaxation order for the Generalized Moment Problem
order = 8 

# polnomial variables
@polyvar x y

# define set K of which the volume should be computed
K = @set(1-x^2-y^2>=0)
# define superset of K where the moments of the Lebegue measure can be computed
B = @set(1-x^2>=0 && 1-y^2>=0)

# define measures on K and B
gmp = GMPModel()
@measure gmp μ [x,y] support=K
@measure gmp ν [x,y] support=B


# We want to maximize the mass of mu
@objective gmp Max Mom(μ,1)

# The constraint is mu + nu = Lebesgue on B
# The normalized moments of Lebesgue measure on B
leb_mom(i,j) = ((1-(-1)^(i+1))/(i+1))*((1-(-1)^(j+1))/(j+1))/4

# DynamicPolynomials.jl provides the possibility to define a monomial vector.
mons = monomials([x,y],0:2*order)
# The monomial vector is not of a polynomial type, we need to convert it first.
pons = polynomial.(mons)

# The fiels mons.Z provides the exponents of the monomials of mons (and pons)
exponents = mons.Z
# The constraints are added to the gmp model. Note that we give a name to the constraints 
# in order to be able to refer to them later on

@constraint gmp cons[i=1:length(exponents)] Mom(pons[i],μ)+Mom(pons[i],ν) == leb_mom(exponents[i]...) 

# In order to relax the problem, we need to specify the relaxation order and a solver factory:
relax!(gmp,order,with_optimizer(Mosek.Optimizer))


# The optimal value is an approximation to the volume of K (which is π in this case)
println("Volume of approximation: $(objective_value(gmp)*4)")

# The polynomial over approximation of the indicator function can be defined
# from the dual solution as follows:
poly = sum(dual_value(gmp,cons[i])*pons[i] for i=1:length(mons))

# Plot 
xx = yy = range(-1, stop = 1, length = 100)
f(xx,yy) = poly(x=>xx,y=>yy)
plot(xx, yy, f, st= :surface)

