#created: TW 2019/02/22
include("../functions.jl")
include("shcl_01.jl")
#include("shcl.jl")

using QuadGK

PDE = Dict()

PDE[:bounds] = Dict()
PDE[:bounds][:t_max]= 1
PDE[:bounds][:x_min]=-1/2
PDE[:bounds][:x_max]= 1/2
PDE[:bounds][:y_min]=0
PDE[:bounds][:y_max]=1

PDE[:flux] = y-> 0.25*y^2
PDE[:initial] = x-> if x<0 return 1 else return 0 end

PDE[:boundaries] = Dict()
#PDE[:boundaries][:periodic] = true
PDE[:boundaries][:left] = t-> 1 
#PDE[:boundaries][:right] = t-> 0

SOL = shcl_01(PDE,2)

for i = 1: length(SOL[:moments][:mu][:moments])
	println((SOL[:moments][:mu][:basis][i],SOL[:moments][:mu][:moments][i]))
end






