include("functions.jl")
using Plots

order = 8

m_over = MomentModel(solver = MosekSolver())
m_linfty = MomentModel(solver = MosekSolver())

# polnomial variable
@polyvar x y

#q = x^2*(x-1/10*y)*(y^2+x)
#q = x*y

q =(x+0.5)*y + (0.5-x)*y^2
K = @set(1-x^2>=0 && 1-y^2>=0)

#m_over

add_measure!(m_over, :mu)
add_support!(m_over, :mu, K)

h(int) = Polynomial{true,Float64}(x^(int[1]))
leb_mom(int) = ((1-(-1)^(int[1]+1))/(int[1]+1))/2
add_mfun_constraint!(m_over, [(h,:mu)], :eq, leb_mom,1)
add_objective!(m_over, :Max, [(Polynomial{true,Float64}(q),:mu)])
mm_solve(m_over,order)

p_over = sum(m_over[:sosmodel].colVal[i+1]*h(i) for i in 0:2*order)

# linfty

add_measure!(m_linfty, :mu)
add_support!(m_linfty, :mu, K)
add_measure!(m_linfty, :nu)
add_support!(m_linfty, :nu, K)

hp(int) = Polynomial{true,Float64}(x^(int[1]))
hm(int) = Polynomial{true,Float64}(-x^(int[1]))
zvec(int) = 0

add_mfun_constraint!(m_linfty, [(hp,:mu),(hm,:nu)], :eq, zvec,1)
add_mconstraint!(m_linfty,[([Polynomial{true,Float64}(1)],:mu),([Polynomial{true,Float64}(1)],:nu)],:eq, [1])
add_objective!(m_linfty, :Max, [(Polynomial{true,Float64}(q),:mu),(Polynomial{true,Float64}(-q),:nu)])
mm_solve(m_linfty,order)
p_linfty = sum(m_linfty[:sosmodel].colVal[i+1]*hp(i) for i in 0:2*order)


X = collect(-1:0.01:1)
OVER = zeros(size(X))
LINF = zeros(size(X))
PROJ = zeros(size(X))
PROJmin = ones(size(X))
Y = collect(-1:0.01:1)
for i = 1:length(X)
    OVER[i] = p_over(x=>X[i])
    LINF[i] = p_linfty(x=>X[i])
    for j = 1:length(Y)
        PROJ[i] = max(q(x=>X[i],y=>Y[j]),PROJ[i])
        PROJmin[i] = min(q(x=>X[i],y=>Y[j]),PROJmin[i])
    end
end

#=
pyplot()
plot(X,PROJ,color = :black,title=order) # projection
plot!(X,OVER, color=:blue)  # over approximation
plot!(X,LINF,color = :green) # linfty approximation
plot!(X,PROJmin,color = :black, style =:dot) # linfty approximation
=#
