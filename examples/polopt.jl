"""
 Solve a global optimization problem with MomentOpt.jl 
 For theortic background see : https://epubs.siam.org/doi/10.1137/S1052623400366802

 Let K be a basic semialgebraic compact set and f a polynomial function. We are interested in  the global optimization problem
 
 min f(x)
     x∈ K.

An equivalent formulation of this problem as Generalized Moment Problem is given by

min <μ,f>
    <μ,1> = 1
     μ∈ M_+(K)

The optimal measure in the Generalized Moment Problem above will be atomic and each atom will be a global minimizer of f on K. 
"""

# Load MomentOpt
using MomentOpt

# Load an implementation of MultivariatePolynomials. Here we use DynamicPolynomials.
using TypedPolynomials

# In order to solve the optimization problem we also need to load an SDP solver.
using MosekTools

# Define polynomial variables
@polyvar x y

# Polynomial to optimize 
f = x^4*y^2 + x^2*y^4 -3x^2*y^2 + 1 

# Define semi algebraic support for the measure 
K = @set(1-x^2 >=0 && 1-y^2 >=0)

gmp = GMPModel()
# Add a variable measure to the model
@variable gmp μ Meas([x,y], support = K)

# Define the objective 
@objective gmp Min Mom(f, μ)
 
# Constrain μ to be a probablity measure
@constraint gmp Mom(1, μ) == 1

# Set the optimizer to be our SDP-solver
set_optimizer(gmp, Mosek.Optimizer)

# The optimize! call generates a relaxation an solves it.
optimize!(gmp)

# By default the relaxation degree is the maximal degree of the registered data. In this case the objective function Mom(f, μ) has degree 6.
println("Relaxation degree: 6")
println("Lower bound: $(objective_value(gmp))")

# We try to extract atoms from the relaxed moment sequence of μ
opt = atomic( μ, tol = 1e-03)
println()

# As we could not extract atoms from the solution, we increase the relaxation degree

set_approximation_degree(gmp, 8)
optimize!(gmp)

println("Relaxation degree: 8")
println("Lower bound: $(objective_value(gmp))")
opt = atomic( μ, tol = 1e-03)

# This time the atom extraction succeeds, which proves optimality of the moment relaxation. 
