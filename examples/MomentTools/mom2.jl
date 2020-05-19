using MomentOpt
using DynamicPolynomials
using MosekTools;optimizer = Mosek.Optimizer

@polyvar x1 x2
e1 = x1^2-2
e2 = (x2^2-3)*(x1*x2-2)

p1 = x1
p2 = 2-x2

m = GMPModel()
@variable m μ Meas([x1, x2], support = @set(e1 ==0 && e2 == 0 && p1 >=0 &&  p2 >= 0))
@constraint m Mom(1, μ) == 1
@objective m Max Mom(x1, μ)

set_approximation_degree(m, 8)
set_optimizer(m, optimizer)
optimize!(m)

v = objective_value(m)
M = MomentOpt.approximation_model(m)

Xi = atomic(μ)

if Xi isa Nothing
    println("Using MultivariateSeries for extraction")
    using MultivariateSeries
    mu = measure(μ)
    y = moment_value.(moments(mu))
    L = monomials(mu)
    _, Xi = decompose(sum(series(yi,l) for (l, yi) in zip(L, y)))
end
Xi
