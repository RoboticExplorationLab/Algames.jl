export
	CollisionConstraint,
	BoundaryConstraint,
	ExpCircleConstraint,
	add_collision_avoidance,
	add_leader_constraints,
	add_circle_boundary,
	state_dim,
	evaluate

struct CollisionConstraint{T,P} <: TO.AbstractConstraint{Inequality,State,P}
	n::Int
	radius1::SVector{P,T} # radius of object 1
	radius2::SVector{P,T} # radius of object 2
	x1::Int  # index of x-state of object 1
	y1::Int  # index of y-state of object 1
	x2::Int  # index of x-state of object 2
	y2::Int  # index of y-state of object 2
	CollisionConstraint(n::Int, radius1::SVector{P,T}, radius2::SVector{P,T},
			x1::Int, y1::Int, x2::Int, y2::Int) where {T,P} =
		 new{T,P}(n,radius1,radius2,x1,y1,x2,y2)
end
TO.state_dim(con::CollisionConstraint) = con.n

function TO.evaluate(con::CollisionConstraint{T,P}, x::SVector) where {T,P}
	x1 = con.x1; x2 = con.x2
	y1 = con.y1; y2 = con.y2
	r = con.radius1 + con.radius2
	-(x[x1] .- x[x2]).^2 - (x[y1] .- x[y2]).^2 .+ r.^2
end


struct BoundaryConstraint{T,P} <: TO.AbstractConstraint{Inequality,State,1}
	n::Int
	p1::SVector{P,T} # initial point of the boundary
	p2::SVector{P,T} # final point of the boundary
	v::SVector{P,T} # vector orthogonal to (p2 - p1) and indicating the forbiden halfspace
	#              p1 *
	#                 |///
	#  NO constraint  |///         constraint
	#  violation      |------> v   violation
	#                 |///
	#                 |///
	#              p2 *
	x::Int  # index of x-state
	y::Int  # index of y-state
	BoundaryConstraint(n::Int, p1::SVector{P,T}, p2::SVector{P,T},
			v::SVector{P,T}, x::Int, y::Int) where {T,P} =
		 	new{T,P}(n,p1,p2,v,x,y)
end
TO.state_dim(con::BoundaryConstraint) = con.n

function TO.evaluate(con::BoundaryConstraint{T,P}, x::SVector) where {T,P}
	x_ = @SVector [x[con.x], x[con.y]]
	if ((x_-con.p1)'*(con.p2-con.p1) > 0) && ((x_-con.p2)'*(con.p1-con.p2) > 0)
		return @SVector [(x_-con.p1)'*con.v]
	else
		return @SVector [0.0]
	end
end


struct ExpCircleConstraint{T,P} <: TO.AbstractConstraint{Inequality,State,P}
	n::Int
	x::SVector{P,T}
	y::SVector{P,T}
	radius::SVector{P,T}
	xi::Int  # index of x-state
	yi::Int  # index of y-state
	ExpCircleConstraint(n::Int, xc::SVector{P,T}, yc::SVector{P,T}, radius::SVector{P,T},
			xi=1, yi=2) where {T,P} =
		 new{T,P}(n,xc,yc,radius,xi,yi)
end
TO.state_dim(con::ExpCircleConstraint) = con.n

function TO.evaluate(con::ExpCircleConstraint{T,P}, x::SVector) where {T,P}
	xc = con.x; xi = con.xi
	yc = con.y; yi = con.yi
	r = con.radius
	exp.(-(x[xi] .- xc).^2 - (x[yi] .- yc).^2 + r.^2) .- 1.0
	# c = -(x[xi] .- xc).^2 - (x[yi] .- yc).^2 + r.^2
	# sign(c[1])*sqrt.(abs.(c))
end


function add_collision_avoidance(conSet::ConstraintSet,
    actors_radii::Vector{T}, px::Vector{Vector{Int}}, p::Int, con_inds::UnitRange) where T
	for i = 1:p
        for j = 1:i-1
            radiusi = @SVector fill(actors_radii[i], 1)
            radiusj = @SVector fill(actors_radii[j], 1)
            col_con = CollisionConstraint(conSet.n, radiusi, radiusj,
                px[i][1], px[i][2], px[j][1], px[j][2])
            add_constraint!(conSet, col_con, con_inds)
        end
    end
    return nothing
end

function add_leader_constraints(conSet::ConstraintSet, leader_states::Vector{SVector{n1,T}},
    actors_radii::Vector{T}, px::Vector{Vector{Int}}) where {n1,T}
    for k = 1:N
        radius_leader = actors_radii[1]
        radius_follower = actors_radii[2]
        con = CircleConstraint(conSet.n,
            SVector{1}([leader_states[k][1]]),
            SVector{1}([leader_states[k][2]]),
            SVector{1}([radius_leader+radius_follower]),
            px[1]...)
        add_constraint!(conSet, con, k:k)
    end
    return nothing
end

function add_circle_boundary!(conSet::ConstraintSet, inds::Array{Int,1},
	x1::Array{T,1}, x2::Array{T,1}, radius::T, M::Int, n::Int, N::Int) where T
	x = SVector{M}([x1[1] + (x2[1] - x1[1])*(j-1)/(M-1) for j=1:M])
	y = SVector{M}([x1[2] + (x2[2] - x1[2])*(j-1)/(M-1) for j=1:M])
	radii = @SVector ones(M)
	radii *= radius
	con = CircleConstraint(n, x, y, radii, inds[1], inds[2])
	add_constraint!(conSet, con, 1:N)
	return nothing
end
