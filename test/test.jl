include("GMP.jl")

using Test

using .GMP
using SemialgebraicSets
using DynamicPolynomials
using SumOfSquares
using CSDP

@polyvar x y
K = @set 1-x^2-y^2>=0
B = @set 1-x^2 >=0


μ = Measure("μ",[x,y],K,SOSCone())
gmp = GMPModel()
add_measure!(gmp, μ)	

ν = Measure("ν",[x,y],support = B)
@measure gmp ν
@measure gmp λ [x]



m1 = Mom(μ,x+y)
m2 = Mom(μ,3)
m3 = Mom(λ, 1.5*x^2)

println(typeof([Mom(μ,1),Mom(μ,1.5)]))
println(typeof([m1,m2]))


n1 = Mom(ν,x+y)
n2 = Mom(ν,1)

me1 = n1+n2 

println(me1)

me2 = Mom(ν,x)+(1.5*Mom(μ,x*y))
println(me2)

@test Mom(y,λ)==nothing

mc1 = MomCons(me1,:eq,0)
mcv1 = MomCons([me1,me1],[:leq,:geq],[1,-1])
mcv2 = MomCons([me1,me2],:eq,0)
mcv3 = MomCons([me1,me2],:eq,[0,1.1])

mass = CMom(μ,1)

@mconstraint(gmp, MomCon(mass,:eq, 1))
@mconstraint(gmp,mc1)
@mconstraints(gmp,mcv1)

@mobjective gmp MomObj(:Max, me1)

@mconstraint gmp MomCon(λ,1,:eq,1)

relax!(gmp,2,with_optimizer(CSDP.Optimizer))

