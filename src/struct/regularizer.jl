################################################################################
# Regularizer
################################################################################

@with_kw mutable struct Regularizer{T}

	# Primal regularization
	"Jacobian regularization for primal variables."
	x::T=1e-3
	u::T=1e-3

	# Dual regularization
	"Jacobian regularization for dual variables."
	λ::T=1e-3
end

import Base.map!
function map!(reg::Regularizer, f)
	for name in fieldnames(Regularizer)
		setfield!(reg, name, f(getfield(reg, name)))
	end
	return nothing
end

function set!(reg::Regularizer, val)
	f(x) = val
	map!(reg, f)
	return nothing
end

function mult!(reg::Regularizer, val)
	f(x) = x*val
	map!(reg, f)
	return nothing
end
