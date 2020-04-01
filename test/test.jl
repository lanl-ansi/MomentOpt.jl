using Test
using DynamicPolynomials

@polyvar x

using Revise
using MomentOpt

m = GMPModel()
@test sprint(show, m) isa String # not important what exactly is printed

@variable m μ Meas([x])
@variable m ν Meas([x]; support = @set x*(1-x) ≥ 0)
@test variables(μ) == variables(ν)
@test !(support(μ) == support(ν))

me = μ + ν

@test sprint(show, me) isa String # not important what exactly is printed
@test MomentOpt.vref_type(me) == MomentOpt.AbstractGMPMeasure
@test variables(me) == [x]

@variable m f Cont([x])
@test_throws AssertionError μ + f

@variable m g Cont([x]; domain = @set x*(1-x) ≥ 0)
@test variables(f) == variables(g)
@test !(domain(f) == domain(g))

fe = f + g

@polyvar y

@variable m ϕ Meas([x,y])
@test_throws AssertionError μ + ϕ

p1 = x^2 + 4.0*x - 1
m1 = Mom(p1, μ)
@test sprint(show, m1) isa String # not important what exactly is printed
m2 = Mom(p1, me)
p2 = x^2 + 4.0*y - 1
@test_throws AssertionError Mom(p2, μ)
m3 = Mom(p2, ϕ)

@test Mom.([p1,p2], [me, ϕ]) isa Vector{<:MomentOpt.MomentExpr}

@constraint m mc2 m2 == 0



