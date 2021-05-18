############################################################################################
#                           CYLINDER CONSTRAINT          								   #
############################################################################################
"""
	CylinderConstraint{P,T}

Constraint of the form

	   NO constraint violation
-----------------------------------------
				NO constraint
				violation
				       r
	                 ---->
                 ////|////  ^
                 ////|////  |
  NO constraint  /// v ///  |     NO constraint
  violation      ////^////  |l    violation
                 ////|con.  |
                 ////|vio.  |
                 / p * ///  |

-----------------------------------------
	   NO constraint violation
The p is the origin, v is the direction, l is the length, r is the radius

# Constructor:
```julia
CylinderConstraint(n, p1::SVector{P}, p2::SVector{P}, p3::SVector{P},
	v::SVector{P},  l::SVector{P}, r::SVector{P}, x=1, y=2, z=3)
```
"""
struct CylinderConstraint{P,T} <: TrajectoryOptimization.StateConstraint
	n::Int
	p1::SVector{P,T}
	p2::SVector{P,T}
	p3::SVector{P,T}
	v::SVector{P,Symbol}
	l::SVector{P,T}
	r::SVector{P,T}
	x::Int  # index of x-state
	y::Int  # index of y-state
	z::Int  # index of z-state
	function CylinderConstraint{P,T}(n::Int,
		p1::AbstractVector, p2::AbstractVector, p3::AbstractVector,
		v::AbstractVector, l::AbstractVector, r::AbstractVector,
		x=1, y=2, z=3) where {P,T}
    	@assert length(p1) == length(p2) == length(p3) == length(v) == length(l) == length(r) == P "Lengths of p1, p2, p3, v, l, p must be equal. Got lengths ($(length(x1)), $(length(y1)), $(length(x2)), $(length(y2)), $(length(xv)), $(length(yv)))"
        new{P,T}(n, p1, p2, p3, v, l, r, x, y, z)
    end
end

function CylinderConstraint(n::Int,
		p1::AbstractVector, p2::AbstractVector, p3::AbstractVector,
		v::AbstractVector, l::AbstractVector, r::AbstractVector,
		x=1, y=2, z=3)
    T = promote_type(
		eltype(p1), eltype(p2), eltype(p3),
		eltype(l), eltype(r))
    P = length(p1)
    CylinderConstraint{P,T}(n, p1, p2, p3, v, l, r, x, y, z)
end

TrajectoryOptimization.state_dim(con::CylinderConstraint) = con.n

function TrajectoryOptimization.evaluate(con::CylinderConstraint, X::StaticVector)
	p1 = con.p1
	p2 = con.p2
	p3 = con.p3
	v = con.v
	l = con.l
	r = con.r
	x = X[con.x]
	y = X[con.y]
	z = X[con.z]

	# t0 = x - p
	t0_1 = x .- p1
	t0_2 = y .- p2
	t0_3 = z .- p3

	valid =
		((v .== :x) .& (t0_1 .> 0.0) .& (t0_1 .< l)) .|
		((v .== :y) .& (t0_2 .> 0.0) .& (t0_2 .< l)) .|
		((v .== :z) .& (t0_3 .> 0.0) .& (t0_3 .< l))

	out = r.^2 - t0_1.^2 - t0_2.^2 - t0_3.^2 +
		(v .== :x) .* t0_1.^2 +
		(v .== :y) .* t0_2.^2 +
		(v .== :z) .* t0_3.^2
	return out .* valid
end

function TrajectoryOptimization.jacobian!(∇c, con::CylinderConstraint{P}, X::SVector) where P
	p1 = con.p1
	p2 = con.p2
	p3 = con.p3
	v = con.v
	l = con.l
	r = con.r
	x = X[con.x]
	y = X[con.y]
	z = X[con.z]

	# t0 = x - p
	t0_1 = x .- p1
	t0_2 = y .- p2
	t0_3 = z .- p3

	valid =
		((v .== :x) .& (t0_1 .> 0.0) .& (t0_1 .< l)) .|
		((v .== :y) .& (t0_2 .> 0.0) .& (t0_2 .< l)) .|
		((v .== :z) .& (t0_3 .> 0.0) .& (t0_3 .< l))

	out = - r.^2 + t0_1.^2 + t0_2.^2 + t0_3.^2 -
		(v .== :x) .* t0_1.^2 -
		(v .== :y) .* t0_2.^2 -
		(v .== :z) .* t0_3.^2

	for i = 1:P
		∇c[i,con.x] = - valid[i] * 2*t0_1[i] * !(v[i] == :x)
		∇c[i,con.y] = - valid[i] * 2*t0_2[i] * !(v[i] == :y)
		∇c[i,con.z] = - valid[i] * 2*t0_3[i] * !(v[i] == :z)
	end
	return false
end

@inline Base.length(::CylinderConstraint{P}) where P = P
@inline TrajectoryOptimization.sense(::CylinderConstraint) = Inequality()

function TrajectoryOptimization.change_dimension(con::CylinderConstraint, n::Int, m::Int, ix=1:n, iu=1:m)
	CylinderConstraint(n,
		con.p1, con.p2, con.p3,
		con.v1, con.v2, con.v3,
		con.l, con.r,
		ix[con.x], ix[con.y], ix[con.z])
end



# # t0 = x - p
# t0_1 = x .- p1
# t0_2 = x .- p1
# t0_3 = x .- p1
# # t1 = t0 x v
# t1_1 = t0_2 .* v3 - t0_3 .* b2
# t1_2 = t0_3 .* v1 - t0_1 .* b3
# t1_3 = t0_1 .* v2 - t0_2 .* b1
# # t2 = norm(t1)
# t2 = sqrt.(t1_1.^2 + t1_2.^2 + t1_3.^2)
# far = t2 .> r
#
# out = (r .- t2).^2
# return out .* bottom .* top .* far



#
# tmp    = (x .- p1) .* v1 + (x .- p2) .* v2 + (x .- p3) .* v3
# top    = tmp .> l
# bottom = tmp .< 0.0
#
# # # t0 = x - p
# # t0_1 = x .- p1
# # t0_2 = y .- p2
# # t0_3 = z .- p3
# # # t1 = norm(t0)
# # t1 = sqrt.(t0_1.^2 + t0_2.^2 + t0_3.^2)
# # # t2 = t0'v
# # t2 = tmp
# # # t3 = t0 - t0'v t0/norm(t0)
# # t3_1 = t0_1 - t2 .* t0_1 ./ t1
# # t3_2 = t0_1 - t2 .* t0_2 ./ t1
# # t3_3 = t0_1 - t2 .* t0_3 ./ t1
# # # t4 = t3 x v
# # t4_1 = t3_2 .* v3 - t3_3 .* b2
# # t4_2 = t3_3 .* v1 - t3_1 .* b3
# # t4_3 = t3_1 .* v2 - t3_2 .* b1
# # # t5 = norm(t3)
# # t5 = sqrt.(t3_1.^2 + t3_2.^2 + t3_3.^2)
#
# far = t5 .> r
# out = (r .- t5).^2
