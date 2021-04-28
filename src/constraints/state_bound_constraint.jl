
############################################################################################
# 							STATE BOUND CONSTRAINTS 										   #
############################################################################################
"""
	StateBoundConstraint{P,N,T}

Linear bound constraint on states
# Constructors
```julia
StateBoundConstraint(n; x_min, x_max)
```
Any of the bounds can be ±∞. The bound can also be specifed as a single scalar, which applies the bound to all states.
"""
struct StateBoundConstraint{P,N,T} <: TrajectoryOptimization.StateConstraint
	n::Int
	x_max::SVector{N,T}
	x_min::SVector{N,T}
	i_max::Vector{Int}
	i_min::Vector{Int}
	inds::SVector{P,Int}
end

Base.copy(bnd::StateBoundConstraint{P,n,T}) where {P,n,T} =
	StateBoundConstraint(bnd.n, bnd.x_max, bnd.x_min,
		copy(bnd.i_max), copy(bnd.i_min), bnd.inds)

function StateBoundConstraint(n; x_max=Inf*(@SVector ones(n)), x_min=-Inf*(@SVector ones(n)))

	# Check and convert bounds
	x_max, x_min = checkBounds(n, x_max, x_min)

	# Concatenate bounds
	b = [-x_max; x_min]
	inds = findall(isfinite, b)
	inds = SVector{length(inds)}(inds)

	# Get linear indices of 1s of Jacobian
	a_max = findall(isfinite, x_max)
	a_min = findall(isfinite, x_min)
	u = length(a_max)
	l = length(a_min)
	carts_u = [CartesianIndex(i,   j) for (i,j) in enumerate(a_max)]
	carts_l = [CartesianIndex(i+u, j) for (i,j) in enumerate(a_min)]
	∇c = zeros(u+l, n)
	linds_u = LinearIndices(zeros(u+l,n))[carts_u]
	linds_l = LinearIndices(zeros(u+l,n))[carts_l]

	StateBoundConstraint(n, x_max, x_min, linds_u, linds_l, inds)
end

function TrajectoryOptimization.con_label(con::StateBoundConstraint, ind::Int)
	i = con.inds[ind]
	n = state_dim(con)
	if 1 <= i <= n
		return "x max $i"
	elseif n < i <= 2n
		j = i - n
		return "x min $j"
	else
		throw(BoundsError())
	end
end

@inline TrajectoryOptimization.state_dim(con::StateBoundConstraint) = con.n
@inline TrajectoryOptimization.is_bound(::StateBoundConstraint) = true
@inline TrajectoryOptimization.lower_bound(bnd::StateBoundConstraint) = bnd.x_min
@inline TrajectoryOptimization.upper_bound(bnd::StateBoundConstraint) = bnd.x_max
@inline TrajectoryOptimization.sense(::StateBoundConstraint) = Inequality()
@inline Base.length(con::StateBoundConstraint) = length(con.i_max) + length(con.i_min)

function TrajectoryOptimization.primal_bounds!(zL, zU, bnd::StateBoundConstraint)
	for i = 1:length(zL)
		zL[i] = max(bnd.x_min[i], zL[i])
		zU[i] = min(bnd.x_max[i], zU[i])
	end
	return true
end

function TrajectoryOptimization.evaluate(bnd::StateBoundConstraint, z::AbstractKnotPoint)
	[(state(z) - bnd.x_max); (bnd.x_min - state(z))][bnd.inds]
end

function TrajectoryOptimization.jacobian!(∇c, bnd::StateBoundConstraint{U,L}, z::AbstractKnotPoint) where {U,L}
	for i in bnd.i_max
		∇c[i]  = 1
	end
	for i in bnd.i_min
		∇c[i] = -1
	end
	return true
end

TrajectoryOptimization.∇jacobian!(G, con::StateBoundConstraint, z::AbstractKnotPoint, λ::AbstractVector) = true # zeros

function TrajectoryOptimization.change_dimension(con::StateBoundConstraint, n::Int, ix=1:n)
	n0 = con.n
	x_max = fill(Inf,n)
	x_min = fill(Inf,n)
	x_max[ix] = con.x_max[1:n0]
	x_min[ix] = con.x_min[1:n0]
	StateBoundConstraint(n, x_max=x_max, x_min=x_min)
end
