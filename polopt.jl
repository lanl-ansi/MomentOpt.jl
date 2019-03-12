#= solve a polynomial optimization problem with JuliaMoments
test two polyvar set up =#

using SCS
include("functions.jl")

# define relaxation order
d = 2

# create moment model and specify solver
mm = MomentModel(solver = SCSSolver())

# define polnomial variable
@polyvar x
@polyvar y

# define measure variable
add_measure!(mm, :mu)

# define semi algebraic support for the measure
K = @set(1-x^2-y^2>=0)
add_support!(mm, :mu, K)

# define the objective funtion
add_objective!(mm, :Min, [(x^2-x*y+y^2,:mu)])

# constrain :mu to be a probablity measure
add_mconstraint!(mm, [([convert(Polynomial{true,Float64},1)],:mu)], :eq, [1])

# solve the moment model
mm_solve(mm,d)

# inspect moments
mm[:moments]

