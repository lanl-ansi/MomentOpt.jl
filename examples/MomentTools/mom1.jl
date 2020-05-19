using MomentOpt
using DynamicPolynomials
using MosekTools
optimizer = Mosek.Optimizer

X = @polyvar x y

p1 = x*y-1
p2 = x^2-2
q  = 2*x-x^2

m = GMPModel()
@variable m μ Meas([x, y], support = @set(p1 == 0 && p2 == 0 && q >= 0))
@constraint m Mom(1, μ) == 1
@objective m Min Mom(y, μ)

set_approximation_degree(m, 6)
set_optimizer(m, optimizer)
optimize!(m)


mu = measure(μ)

v = objective_value(m)
y = moment_value.(moments(mu))
L = monomials(mu)

Xi = atomic(μ)

if Xi isa Nothing
    println("Using MultivariateSeries for extraction")
    using MultivariateSeries
    _, Xi = decompose(sum(series(yi,l) for (l,yi) in zip(L, y)))
end
Xi
