############################################################################################
#                              WALL CONSTRAINT          								   #
############################################################################################
"""
	WallConstraint{P,T}

Constraint of the form

	   NO constraint violation
-----------------------------------------
              p1 * ///
                 |////
  NO constraint  |////        constraint
  violation      |------> v   violation
                 |////
                 |////
              p2 * ///
-----------------------------------------
	   NO constraint violation

`` (x - p1)^T v \\leq 0 ``
where ``x = [x[x], x[y]]``, ``p1 = [x1,y1]``, ``p2 = [x2,y2]`` ``v = [xv,yv]``.

# Constructor:
```julia
WallConstraint(n, x1::SVector{P}, y1::SVector{P}, x2::SVector{P}, y2::SVector{P}, xv::SVector{P}, yv::SVector{P}, x=1, y=2)
```
"""
struct WallConstraint{P,T} <: TrajectoryOptimization.StateConstraint
	n::Int
	x1::SVector{P,T}
	y1::SVector{P,T}
	x2::SVector{P,T}
	y2::SVector{P,T}
	xv::SVector{P,T}
	yv::SVector{P,T}
	x::Int  # index of x-state
	y::Int  # index of y-state
	function WallConstraint{P,T}(n::Int, x1::AbstractVector, y1::AbstractVector,
		x2::AbstractVector, y2::AbstractVector, xv::AbstractVector, yv::AbstractVector,
		x=1, y=2) where {P,T}
    	@assert length(x1) == length(y1) == length(x2) == length(y2) == length(xv) == length(yv) == P "Lengths of x1, y1, x2, y2, xv, yv must be equal. Got lengths ($(length(x1)), $(length(y1)), $(length(x2)), $(length(y2)), $(length(xv)), $(length(yv)))"
        new{P,T}(n, x1, y1, x2, y2, xv, yv, x, y)
    end
end
function WallConstraint(n::Int, x1::AbstractVector, y1::AbstractVector,
	x2::AbstractVector, y2::AbstractVector, xv::AbstractVector, yv::AbstractVector,
		x=1, y=2)
    T = promote_type(eltype(x1), eltype(y1), eltype(x2), eltype(y2), eltype(xv), eltype(yv))
    P = length(x1)
    WallConstraint{P,T}(n, x1, y1, x2, y2, xv, yv, x, y)
end
TrajectoryOptimization.state_dim(con::WallConstraint) = con.n

function TrajectoryOptimization.evaluate(con::WallConstraint, X::StaticVector)
	x1 = con.x1
	y1 = con.y1
	x2 = con.x2
	y2 = con.y2
	xv = con.xv
	yv = con.yv
	x = X[con.x]
	y = X[con.y]

	left  = (x .- x1) .* (x2 .- x1) + (y .- y1) .* (y2 .- y1) .> 0
	right = (x .- x2) .* (x1 .- x2) + (y .- y2) .* (y1 .- y2) .> 0
	out = (x .- x1) .* xv + (y .- y1) .* yv
	return out .* left .* right
end

function TrajectoryOptimization.jacobian!(∇c, con::WallConstraint{P}, X::SVector) where P
	x1 = con.x1
	y1 = con.y1
	x2 = con.x2
	y2 = con.y2
	xv = con.xv
	yv = con.yv
	x = X[con.x]
	y = X[con.y]

	left  = (x .- x1) .* (x2 .- x1) + (y .- y1) .* (y2 .- y1) .> 0
	right = (x .- x2) .* (x1 .- x2) + (y .- y2) .* (y1 .- y2) .> 0
	for i = 1:P
		∇c[i,con.x] = left[i] * right[i] * xv[i]
		∇c[i,con.y] = left[i] * right[i] * yv[i]
	end
	return false
end

@inline Base.length(::WallConstraint{P}) where P = P
@inline TrajectoryOptimization.sense(::WallConstraint) = Inequality()

function TrajectoryOptimization.change_dimension(con::WallConstraint, n::Int, m::Int, ix=1:n, iu=1:m)
	WallConstraint(n, con.x1, con.y1, con.x2, con.y2, con.xv, con.yv, ix[con.x], ix[con.y])
end
