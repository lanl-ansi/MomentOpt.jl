#= robust approximation
=#

order = 5
using Plots
include("functions.jl")
mm = MomentModel(solver = MosekSolver())

# polnomial variable
@polyvar x y

add_measure!(mm, :mu)

K = @set(1-x^2>=0&&1-y^2>=0)

add_support!(mm, :mu, K)
# q = 0.25-x^4 + 1/2*y*x^2 +1/2*x*y^2 -y^4
# q = 0.25-x^2 -1/2*y^2 + x*y
q = x^2*(x-1/10*y)*(y^2+x)

add_objective!(mm, :Max, [(Polynomial{true,Float64}(q),:mu)])


# maj constraint
h(int) = Polynomial{true,Float64}(x^(int[1]))
# normalized moments of Lebesgue measure on B
leb_mom(int) = ((1-(-1)^(int[1]+1))/(int[1]+1))/2

add_mfun_constraint!(mm, [(h,:mu)], :eq, leb_mom,1)

mm_solve(mm,order)

p = sum(mm[:sosmodel].colVal[i]*h(i-1) for i in 1:2*order+1)
 
#=
# plotting

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
wireframe(XX,YY,ZZ, color=:blue)
wireframe!(XX,YY,ZZZ, color= :black)
=#
