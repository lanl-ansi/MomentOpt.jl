export Measure, variables, support, certificate

#TODO: problem with type inference (parametrized measures instead Measure)

mutable struct Measure{V<:MP.AbstractVariable,SAS<:AbstractBasicSemialgebraicSet,C<:PJ.PolynomialSet} #,PB<:PJ.AbstractPolynomialBasis}
	name:: String	# name
	vars:: Vector{V}# Vector of variables associated to the measure
	supp:: SAS	# Basic semialgebraic support of the measure
	cert:: C   	# Certificate set like SOSCone(), should depend on a basis, and structure information
#	mons:: PB   	# Basis in which the monomials are expressed
#TODO: Find solution for the basis in which moments are expressed.
end

Base.broadcastable(mu::Measure) = Ref(mu)

# constructor
function Measure(name::AbstractString, vars::Vector{V}; support = FullSpace(), certificate = SOSCone() ) where V<: MP.AbstractVariable
	return Measure(name,vars,support,certificate)
end

variabletype(mu::Measure{V}) where {V} = V

function MP.variables(mu::Measure)
	return mu.vars
end

function support(mu::Measure)
	return mu.supp
end

function certificate(mu::Measure)
	return mu.cert
end

# printing
function Base.show(io::IO, meas::Measure)
	print(io, meas.name)
end
