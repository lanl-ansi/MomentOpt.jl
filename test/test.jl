using DynamicPolynomials

@polyvar x

using Revise
using MomentOpt

struct Dummy1 <: MomentOpt.AbstractApproximation end

mu = VariableMeasure([x], FullSpace(), Dummy1(), MonomialBasis)
nu = VariableMeasure([x], FullSpace(), Dummy1(), MonomialBasis)

m = GMPModel()

@variable m μ mu
@variable m ν nu

me = μ + ν

vf = VariableContinuous(MomentOpt.MOI.get.(mu, MomentOpt.GenericObjectAttributes)...)
vg = VariableContinuous(MomentOpt.MOI.get.(mu, MomentOpt.GenericObjectAttributes)...)

@variable m f vf
@variable m g vg

fe = f + g

@polyvar y

phi =  VariableMeasure([x, y], FullSpace(), Dummy1(), MonomialBasis)

@variable m ϕ phi

