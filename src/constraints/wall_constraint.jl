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
The ordering of p1 and p2 does not matter.

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



############################################################################################
#                              WALL3D CONSTRAINT          								   #
############################################################################################
"""
	Wall3DConstraint{P,T}

Constraint of the form

	   NO constraint violation
-----------------------------------------
  *  -------- p1 * ///
                 |////
  NO constraint  |////        constraint
  violation      |------> v   violation
                 |////
                 |////
p3 * -------- p2 * ///

                 -
	         -|  NO cons. viol.
      p3 *c.v.|
     -    |   |               -
 - |  cons.|  |           -
   |  viol. | |       -
   |         ||   -
p1 * -------- * p2
   |  c.v. -
   |  -
  -|     NO cons. viol.
-----------------------------------------
	   NO constraint violation
The points p1, p2, p3 form a parallelepiped in 3D space.

`` (x - p1)^T v \\leq 0 ``
where ``x = [x[x], x[y], x[z]]``, ``p1 = [x1,y1,z1]``, ``p2 = [x2,y2,z2]``, ``p3 = [x3,y3,z3]``, ``v = [xv,yv,zv]``.

# Constructor:
```julia
Wall3DConstraint(n, x1::SVector{P}, y1::SVector{P}, z1::SVector{P}, x2::SVector{P}, y2::SVector{P}, z2::SVector{P},  x3::SVector{P}, y3::SVector{P}, z3::SVector{P}, xv::SVector{P}, yv::SVector{P}, x=1, y=2, z=3)
```
"""
struct Wall3DConstraint{P,T} <: TrajectoryOptimization.StateConstraint
	n::Int
	x1::SVector{P,T}
	y1::SVector{P,T}
	z1::SVector{P,T}
	x2::SVector{P,T}
	y2::SVector{P,T}
	z2::SVector{P,T}
	x3::SVector{P,T}
	y3::SVector{P,T}
	z3::SVector{P,T}
	xv::SVector{P,T}
	yv::SVector{P,T}
	zv::SVector{P,T}
	x::Int  # index of x-state
	y::Int  # index of y-state
	z::Int  # index of z-state
	function Wall3DConstraint{P,T}(n::Int,
		x1::AbstractVector, y1::AbstractVector, z1::AbstractVector,
		x2::AbstractVector, y2::AbstractVector, z2::AbstractVector,
		x3::AbstractVector, y3::AbstractVector, z3::AbstractVector,
		xv::AbstractVector, yv::AbstractVector, zv::AbstractVector,
		x=1, y=2, z=3) where {P,T}
    	@assert length(x1) == length(y1) == length(z1) == length(x2) == length(y2) == length(z2) == length(x3) == length(y3) == length(z3) == length(xv) == length(yv) == length(zv) == P "Lengths of x1, y1, x2, y2, xv, yv must be equal. Got lengths ($(length(x1)), $(length(y1)), $(length(x2)), $(length(y2)), $(length(xv)), $(length(yv)))"
        new{P,T}(n, x1, y1, z1, x2, y2, z2, x3, y3, z3, xv, yv, zv, x, y, z)
    end
end
function Wall3DConstraint(n::Int,
		x1::AbstractVector, y1::AbstractVector, z1::AbstractVector,
		x2::AbstractVector, y2::AbstractVector, z2::AbstractVector,
		x3::AbstractVector, y3::AbstractVector, z3::AbstractVector,
		xv::AbstractVector, yv::AbstractVector, zv::AbstractVector,
		x=1, y=2, z=3)
    T = promote_type(
		eltype(x1), eltype(y1), eltype(z1),
		eltype(x2), eltype(y2), eltype(z2),
		eltype(x3), eltype(y3), eltype(z3),
		eltype(xv), eltype(yv), eltype(zv))
    P = length(x1)
    Wall3DConstraint{P,T}(n, x1, y1, z1, x2, y2, z2, x3, y3, z3, xv, yv, zv, x, y, z)
end
TrajectoryOptimization.state_dim(con::Wall3DConstraint) = con.n

function TrajectoryOptimization.evaluate(con::Wall3DConstraint, X::StaticVector)
	x1 = con.x1
	y1 = con.y1
	z1 = con.z1
	x2 = con.x2
	y2 = con.y2
	z2 = con.z2
	x3 = con.x3
	y3 = con.y3
	z3 = con.z3
	xv = con.xv
	yv = con.yv
	zv = con.zv
	x = X[con.x]
	y = X[con.y]
	z = X[con.z]

	left    = (x .- x1) .* (x2 .- x1) + (y .- y1) .* (y2 .- y1) + (z .- z1) .* (z2 .- z1) .> 0
	right   = (x .- x2) .* (x1 .- x2) + (y .- y2) .* (y1 .- y2) + (z .- z2) .* (z1 .- z2) .> 0
	bottom  = (x .- x3) .* (x2 .- x3) + (y .- y3) .* (y2 .- y3) + (z .- z3) .* (z2 .- z3) .> 0
	top     = (x .- x2) .* (x3 .- x2) + (y .- y2) .* (y3 .- y2) + (z .- z2) .* (z3 .- z2) .> 0
	out = (x .- x1) .* xv + (y .- y1) .* yv + (z .- z1) .* zv
	return out .* left .* right .* bottom .* top
end

function TrajectoryOptimization.jacobian!(∇c, con::Wall3DConstraint{P}, X::SVector) where P
	x1 = con.x1
	y1 = con.y1
	z1 = con.z1
	x2 = con.x2
	y2 = con.y2
	z2 = con.z2
	x3 = con.x3
	y3 = con.y3
	z3 = con.z3
	xv = con.xv
	yv = con.yv
	zv = con.zv
	x = X[con.x]
	y = X[con.y]
	z = X[con.z]

	left    = (x .- x1) .* (x2 .- x1) + (y .- y1) .* (y2 .- y1) + (z .- z1) .* (z2 .- z1) .> 0
	right   = (x .- x2) .* (x1 .- x2) + (y .- y2) .* (y1 .- y2) + (z .- z2) .* (z1 .- z2) .> 0
	bottom  = (x .- x3) .* (x2 .- x3) + (y .- y3) .* (y2 .- y3) + (z .- z3) .* (z2 .- z3) .> 0
	top     = (x .- x2) .* (x3 .- x2) + (y .- y2) .* (y3 .- y2) + (z .- z2) .* (z3 .- z2) .> 0
	for i = 1:P
		∇c[i,con.x] = left[i] * right[i] * bottom[i] * top[i] * xv[i]
		∇c[i,con.y] = left[i] * right[i] * bottom[i] * top[i] * yv[i]
		∇c[i,con.z] = left[i] * right[i] * bottom[i] * top[i] * zv[i]
	end
	return false
end

@inline Base.length(::Wall3DConstraint{P}) where P = P
@inline TrajectoryOptimization.sense(::Wall3DConstraint) = Inequality()

function TrajectoryOptimization.change_dimension(con::Wall3DConstraint, n::Int, m::Int, ix=1:n, iu=1:m)
	Wall3DConstraint(n,
		con.x1, con.y1, con.z1,
		con.x2, con.y2, con.z2,
		con.x3, con.y3, con.z3,
		con.x4, con.y4, con.z4,
		con.xv, con.yv, con.zv,
		ix[con.x], ix[con.y], ix[con.z])
end
