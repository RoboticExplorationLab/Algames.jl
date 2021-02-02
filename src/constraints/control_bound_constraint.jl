
############################################################################################
# 				     	CONTROL BOUND CONSTRAINTS 										   #
############################################################################################
"""
	ControlBoundConstraint{P,M,T}

Linear bound constraint on controls
# Constructors
```julia
ControlBoundConstraint(m; u_min, u_max)
```
Any of the bounds can be ±∞. The bound can also be specifed as a single scalar, which applies the bound to all controls.
"""
struct ControlBoundConstraint{P,M,T} <: TrajectoryOptimization.ControlConstraint
	m::Int
	u_max::SVector{M,T}
	u_min::SVector{M,T}
	i_max::Vector{Int}
	i_min::Vector{Int}
	inds::SVector{P,Int}
end

Base.copy(bnd::ControlBoundConstraint{P,m,T}) where {P,m,T} =
	ControlBoundConstraint(bnd.m, bnd.u_max, bnd.u_min,
		copy(bnd.i_max), copy(bnd.i_min), bnd.inds)

function ControlBoundConstraint(m; u_max=Inf*(@SVector ones(m)), u_min=-Inf*(@SVector ones(m)))

	# Check and convert bounds
	u_max, u_min = checkBounds(m, u_max, u_min)

	# Concatenate bounds
	b = [-u_max; u_min]
	inds = findall(isfinite, b)
	inds = SVector{length(inds)}(inds)

	# Get linear indices of 1s of Jacobian
	a_max = findall(isfinite, u_max)
	a_min = findall(isfinite, u_min)
	u = length(a_max)
	l = length(a_min)
	carts_u = [CartesianIndex(i,   j) for (i,j) in enumerate(a_max)]
	carts_l = [CartesianIndex(i+u, j) for (i,j) in enumerate(a_min)]
	∇c = zeros(u+l, m)
	linds_u = LinearIndices(zeros(u+l,m))[carts_u]
	linds_l = LinearIndices(zeros(u+l,m))[carts_l]

	ControlBoundConstraint(m, u_max, u_min, linds_u, linds_l, inds)
end

function TrajectoryOptimization.con_label(con::ControlBoundConstraint, ind::Int)
	i = con.inds[ind]
	m = control_dim(con)
	if 1 <= i <= m
		return "u max $i"
	elseif m < i <= 2m
		j = i - m
		return "u min $j"
	else
		throw(BoundsError())
	end
end

function checkBounds(n::Int, u::AbstractVector, l::AbstractVector)
	if all(u .>= l)
		return SVector{n}(u), SVector{n}(l)
	else
		throw(ArgumentError("Upper bounds must be greater than or equal to lower bounds"))
	end
end

checkBounds(n::Int, u::Real, l::Real) =
	checkBounds(n, (@SVector fill(u,n)), (@SVector fill(l,n)))
checkBounds(n::Int, u::AbstractVector, l::Real) = checkBounds(n, u, fill(l,n))
checkBounds(n::Int, u::Real, l::AbstractVector) = checkBounds(n, fill(u,n), l)


@inline TrajectoryOptimization.control_dim(con::ControlBoundConstraint) = con.m
@inline TrajectoryOptimization.is_bound(::ControlBoundConstraint) = true
@inline TrajectoryOptimization.lower_bound(bnd::ControlBoundConstraint) = bnd.u_min
@inline TrajectoryOptimization.upper_bound(bnd::ControlBoundConstraint) = bnd.u_max
@inline TrajectoryOptimization.sense(::ControlBoundConstraint) = Inequality()
@inline Base.length(con::ControlBoundConstraint) = length(con.i_max) + length(con.i_min)

function TrajectoryOptimization.primal_bounds!(zL, zU, bnd::ControlBoundConstraint)
	for i = 1:length(zL)
		zL[i] = max(bnd.u_min[i], zL[i])
		zU[i] = min(bnd.u_max[i], zU[i])
	end
	return true
end

function TrajectoryOptimization.evaluate(bnd::ControlBoundConstraint, z::AbstractKnotPoint)
	[(control(z) - bnd.u_max); (bnd.u_min - control(z))][bnd.inds]
end

function TrajectoryOptimization.jacobian!(∇c, bnd::ControlBoundConstraint{U,L}, z::AbstractKnotPoint) where {U,L}
	for i in bnd.i_max
		∇c[i]  = 1
	end
	for i in bnd.i_min
		∇c[i] = -1
	end
	return true
end

TrajectoryOptimization.∇jacobian!(G, con::ControlBoundConstraint, z::AbstractKnotPoint, λ::AbstractVector) = true # zeros

function TrajectoryOptimization.change_dimension(con::ControlBoundConstraint, m::Int, iu=1:m)
	m0 = con.m
	u_max = fill(Inf,m)
	u_min = fill(Inf,m)
	u_max[iu] = con.u_max[1:m0]
	u_min[iu] = con.u_min[1:m0]
	ControlBoundConstraint(m, u_max=u_max, u_min=u_min)
end
