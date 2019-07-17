"""
 Solve a global optimization problem with MomentOpt.jl 
 For theortic background see : https://epubs.siam.org/doi/10.1137/S1052623400366802

 Let K be a basic semialgebraic compact set and f a polynomial function. We are intersted in 
 the global optimization problem
 
 min f(x)
     x∈K

An equivalent formulation of this problem as Generalized Moment Problem is given by

min <μ,f>
    <μ,1> = 1
     μ∈M(K)

The optimal measure in the Generalized Moment Problem above will be atomic and each atom will be a global minimizer 
of f on K. 
"""

using DynamicPolynomials
using SemialgebraicSets
using SumOfSquares 

using MomentOpt

using MosekTools

# Define polnomial variables
@polyvar x y

# Polynomial to optimize 
f =x^4*y^2 + x^2*y^4 -3x^2*y^2 + 1 

# Define semi algebraic support for the measure 
K = @set(1-x>=0 && x+1>=0 && 1-y>=0 && y+1>=0)


gmp = GMPModel()
# Add measure to the model
@measure gmp μ [x,y] support=K

# Define the objective 
@objective gmp Min Mom(f,μ)
# This notation is a shorthand for 
# @mobjective gmp :Min Mom(f,μ)
 

# Constrain μ to be a probablity measure
@mconstraint gmp MomCon(Mom(1,μ), :eq, 1)


# We solve the relaxation of order 2 with CSDP
relax!(gmp, 2, with_optimizer(Mosek.Optimizer))


println("Relaxation order: $(2)")
println("Objective value: $(objective_value(gmp))")
# We try to extract atoms from the relaxed moment sequence of μ
opt = atomic(gmp,μ)
println()

# As we could not extract atoms from the solution, we increase the relaxation order
relax!(gmp, 3, with_optimizer(Mosek.Optimizer))

println("Relaxation order: $(3)")
println("Objective value: $(objective_value(gmp))")
opt = atomic(gmp,μ)

# This time the atom extraction succeeds, which proves optimality of the moment relaxation. 



