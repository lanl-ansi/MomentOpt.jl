# See `examples/burgers.jl`

using DynamicPolynomials
using SemialgebraicSets

using MomentOpt

function burgers(order)
    f(y) = 1/4 * y^2
    T = [0, 1]
    X = [-1/2, 1/2]
    Y = [0, 1]

    # Define polynomial variables
    @polyvar t x y

    gmp = GMPModel()
    # Add measures

    @measure gmp μ  [t,x,y] support=@set((t-T[1])*(T[2]-t)>=0 && (x-X[1])*(X[2]-x)>=0 && (y-Y[1])*(Y[2]-y)>=0)
    @measure gmp μT [t,x,y] support=@set(t==1 && (x-X[1])*(X[2]-x)>=0 && (y-Y[1])*(Y[2]-y)>=0)
    @measure gmp μR [t,x,y] support=@set((t-T[1])*(T[2]-t)>=0 && x==1/2 && (y-Y[1])*(Y[2]-y)>=0)

    # As T×X×Y is compact we can use the monomials as testfunctions φ. 
    # We define the values for the right hand side of the constraint of the dynamic
    rhs(i,j) = 0^i*(-(-1/2)^(j+1))/(j+1)+ (-1/2)^j/((i+1)*4)

    # DynamicPolynomials.jl provides the possibility to define a monomial vector.
    mons = monomials([t,x],0:2*order-1)
    # The monomial vector is not of a polynomial type, we need to convert it first.
    pons = polynomial.(mons)
    wpde = vec(differentiate(pons,[t]))*y + vec(differentiate(pons,[x]))*f(y) 

    # The fiels mons.Z provides the exponents of the monomials of mons (and pons)
    @constraint gmp Mom.(pons*y,μT)+Mom.(pons*f(y),μR)-Mom.(wpde,μ) .== [rhs(mons.Z[i]...) for i = 1:length(mons.Z)] 

    # Right hand sides for the marginal constraints
    lebt(i) = (T[2]^(i+1)-T[1]^(i+1))/(i+1)
    lebx(j) = (X[2]^(j+1)-X[1]^(j+1))/(j+1)
    lebtx(i,j) = lebt(i)*lebx(j)
    monstx = monomials([t,x], 0:2*order)
    monsx= monomials(x, 0:2*order)
    monst= monomials(t, 0:2*order)

    @constraint gmp Mom.(polynomial.(monstx),μ) .== [lebtx(monstx.Z[i]...) for i =1:length(monstx.Z)]
    @constraint gmp Mom.(polynomial.(monsx),μT) .== [lebx(monsx.Z[i]...) for i =1:length(monsx.Z)]
    @constraint gmp Mom.(polynomial.(monst),μR) .== [lebt(monst.Z[i]...) for i =1:length(monst.Z)]

    # The GMP formulation of the PDE does not have any objective function as the measures are uniquely determined by the
    # problem. However, when relaxing the problem the measures will not be unique anymore which is why we add an objective function.
    # Note that a GMPModel always needs an objective to be relaxed.
    # A common heuristic for an objective function is to minimize the trace of the moment matrix:

    mmons = polynomial.(monomials([t, x, y], 0:order))
    trace = mmons' * mmons

    @objective(gmp, Min, Mom(trace, μ) + Mom(trace, μT) + Mom(trace, μR))

    @time relax!(gmp,order, with_optimizer(() -> MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())))

    return
end
