#= how does q look like =#

order = 6;

using Plots

include("functions.jl")

mm = MomentModel(solver = MosekSolver())

# polnomial variable
@polyvar x

add_measure!(mm, :mu)
add_measure!(mm, :nu)

K = @set(0.25-x^2>=0)
B = @set(1-x^2>=0)

add_support!(mm, :mu, K)
add_support!(mm, :nu, B)

add_objective!(mm, :Max, [(Polynomial{true,Float64}(1),:mu)])


# maj constraint
h(int) = Polynomial{true,Float64}(x^(int[1]))
# normalized moments of Lebesgue measure on B
leb_mom(int) = ((1-(-1)^(int[1]+1))/(int[1]+1))/2

add_mfun_constraint!(mm, [(h,:mu), (h,:nu)], :eq, leb_mom,1)

# stokes constraint
s(int) = Polynomial{true,Float64}((int[1])*x^(int[1]-1)*(0.25-x^2) + x^(int[1])*(-2*x))
# normalized moments of Lebesgue measure on B
null_rhs(int) = 0
add_mfun_constraint!(mm, [(s,:mu)], :eq, null_rhs,1)

mm_solve(mm,order)

mm[:moments]

# plotting
X = -1:0.01:1
PX = zeros(length(X),1)
QX = zeros(length(X),1)
#QK = zeros(length(X),1)
for i = 1:length(X)
    PX[i] = sum(mm[:sosmodel].colVal[j]*h(j-1)(x=>X[i]) for j=1:2*order+1)
    QX[i] = sum(mm[:sosmodel].colVal[j+2*order+1]*s(j-1)(x=>X[i]) for j=1:2*order)
end


# how to choose these variables? find the right indexing
println("Volume is")
println(getobjectivevalue(mm[:sosmodel]))
#=
pyplot()
foo = plot([-1,-0.5,-0.5,0.5,0.5,1],[0,0,1,1,0,0],c=:black,linestyle= :dash,label = :indicator,linewidth=2)
plot!(foo,X,PX,label="p",linewidth=3, c=:red)
plot!(foo,X,QX,label = "q",linewidth=3,c=:blue)
=#
