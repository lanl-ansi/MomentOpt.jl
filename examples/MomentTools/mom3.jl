using MomentOpt
using DynamicPolynomials
using MosekTools
optimizer = Mosek.Optimizer

@polyvar x y
q1 = 1-x^2-y^2
q2 = x^3-y^2

m = GMPModel()
@variable m μ Meas([x, y], support = @set(q1 >=0 &&  q2 >= 0))
@constraint m Mom(1, μ) == 1
@objective m Min Mom(x, μ)

set_approximation_degree(m, 8)
set_optimizer(m, optimizer)
optimize!(m)

v = objective_value(m)
M = MomentOpt.approximation_model(m)

println("Extraction MultivariateMoments")
Xi = atomic(μ)
println(Xi)


println("Extraction MultivariateSeries")
Xi = let
    using MultivariateSeries
    mu = measure(μ)
    y = moment_value.(moments(mu))
    L = monomials(mu)
    _, Xi = decompose(sum(series(yi,l) for (l, yi) in zip(L, y)))
    Xi
end
println(Xi)


"""
MomentTools uses smaller localization matrices. 
MomentTools uses monomials up to degree (degree - maxdegree(p)) to index the localization matrix. The maximal degree in the localization matrix therefore is 2*(degree - maxdegree(p)) + maxdegree(p) = 2*degree- maxdegree(p)

MomentOpt uses monomials up to degree floor(2*degree - maxdegree(p)) and ends up with a maximal degree of 2*degree or 2*degree-1 depending on pairity od maxdegree(p). 

Note that the approximation degree in MomentOpt corresonds to the maximum degree of moments computed, while in MomentTools it corresponds to the maximum degree used for indexing the moment matrix. In the note above we used the terminology of MomentTools. 
"""
