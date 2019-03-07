#= robust approximation
=#

order = 4
using Plots
include("functions.jl")
mm = MomentModel(solver = MosekSolver())

# polnomial variable
@polyvar x
@polyvar y

add_measure!(mm, :mu)
add_measure!(mm, :nu)

K = @set(1-x^2>=0 && 1-y^2>=0)
add_support!(mm, :mu, K)
add_support!(mm, :nu, K)

# maj constraint
hp(int) = Polynomial{true,Float64}(x^(int[1]))
hm(int) = Polynomial{true,Float64}(-x^(int[1]))
zvec(int) = 0

add_mfun_constraint!(mm, [(hp,:mu),(hm,:nu)], :eq, zvec,1)
add_mconstraint!(mm,[([Polynomial{true,Float64}(1)],:mu),([Polynomial{true,Float64}(1)],:nu)],:eq, [1])

#q = 0.25*x^4-x^2 +1/10*x^2*y^2 +1/4*x*3 -1/2*y^2 + 1/4*x*y
q = x^2*(x-1/10*y)*(y^2+x)

add_objective!(mm, :Max, [(Polynomial{true,Float64}(q),:mu),(Polynomial{true,Float64}(-q),:nu)])


mm_solve(mm,order)

mm[:moments]

p = sum(mm[:sosmodel].colVal[i]*x^(i-1) for i in 1:2*order+1)
#= plotting


X = collect(-1:0.1:1)
Y = collect(-1:0.1:1)

n = length(X)
XX = zeros(n,n)
YY = zeros(n,n)
for i = 1:n
    XX[:,i] = X
    YY[i,:] = Y
end


ZZ = zeros(n,n)
ZZZ = zeros(n,n)
for i in 1:n
    for j in 1:n
        ZZ[i,j] = p(x=>X[i])
        ZZZ[i,j] =q(x=>X[i],y => Y[j])
    end
end
pyplot()
wireframe(XX,YY,ZZ)
wireframe!(XX,YY,ZZZ)
=#
