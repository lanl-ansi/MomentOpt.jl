#= solve a polynomial optimization problem with JuliaMoments
test two cosntraint set up
test mfun-constraint =#

order = 6
using SDPA
using Plots

include("functions.jl")

mm = MomentModel(solver = SDPASolver())

# polnomial variable
@polyvar x
@polyvar y

add_measure!(mm, :mu)
add_measure!(mm, :nu)

K = @set(1-x^2-y^2>=0)
B = @set(1-x^2>=0 && 1-y^2>=0)

add_support!(mm, :mu, K)
add_support!(mm, :nu, B)

add_objective!(mm, :Max, [(convert(Polynomial{true,Float64},1),:mu)])


# maj constraint
h(i,j) = convert(Polynomial{true,Float64},x^i*y^j)
# normalized moments of Lebesgue measure on B
leb_mom(i,j) = ((1-(-1)^(i+1))/(i+1))*((1-(-1)^(j+1))/(j+1))/4

add_mfun_constraint!(mm, [(h,:mu), (h,:nu)], :eq, leb_mom,2)

# stokes constraint
sx(i,j) = differentiate(x^i*y^j*(1-x^2-y^2),x)
sy(i,j) = differentiate(x^i*y^j*(1-x^2-y^2),y)

# normalized moments of Lebesgue measure on B
null_rhs(i,j) = 0
add_mfun_constraint!(mm, [(sx,:mu)], :eq, null_rhs,2)
add_mfun_constraint!(mm, [(sy,:mu)], :eq, null_rhs,2)

mm_solve(mm,order)

println("Volume is")
println(4*getobjectivevalue(mm[:sosmodel]))
