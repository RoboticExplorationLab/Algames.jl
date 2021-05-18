################################################################################
# Add Collision Avoidance
################################################################################

function add_collision_avoidance!(game_con::GameConstraintValues, i::Int, j::Int, radius::T) where {T}
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	px = probsize.px

	add_constraint!(game_con.state_conlist[i], CollisionConstraint(n,px[i],px[j],radius), 2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end

function add_collision_avoidance!(game_con::GameConstraintValues, radius::Vector{T}) where {T}
	probsize = game_con.probsize
	p = probsize.p
	@assert p == length(radius)
	for i = 1:p
		for j in setdiff(1:p,i)
			ri = radius[i]
			rj = radius[j]
			add_collision_avoidance!(game_con, i, j, ri+rj)
		end
	end
	return nothing
end

function add_collision_avoidance!(game_con::GameConstraintValues, radius::T) where {T}
	p = game_con.probsize.p
	add_collision_avoidance!(game_con, radius*ones(p))
	return nothing
end

################################################################################
# Add Spherical Collision Avoidance
################################################################################

function add_spherical_collision_avoidance!(game_con::GameConstraintValues, i::Int, j::Int, radius::T) where {T}
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	pz = probsize.pz

	pzi = pz[i][1:3] # assumes that the 3 first dimensions of the state of player i correspond to x,y,z.
	pzj = pz[j][1:3] # assumes that the 3 first dimensions of the state of player i correspond to x,y,z.
	add_constraint!(game_con.state_conlist[i], CollisionConstraint(n,pzi,pzj,radius), 2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end

function add_spherical_collision_avoidance!(game_con::GameConstraintValues, radius::Vector{T}) where {T}
	probsize = game_con.probsize
	p = probsize.p
	@assert p == length(radius)
	for i = 1:p
		for j in setdiff(1:p,i)
			ri = radius[i]
			rj = radius[j]
			add_spherical_collision_avoidance!(game_con, i, j, ri+rj)
		end
	end
	return nothing
end

function add_spherical_collision_avoidance!(game_con::GameConstraintValues, radius::T) where {T}
	p = game_con.probsize.p
	add_spherical_collision_avoidance!(game_con, radius*ones(p))
	return nothing
end

################################################################################
# Add State Bounds
################################################################################

function add_state_bound!(game_con::GameConstraintValues, i::Int, x_max::AbstractVector, x_min::AbstractVector)
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	add_constraint!(game_con.state_conlist[i], StateBoundConstraint(n,x_max=x_max,x_min=x_min), 2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end

################################################################################
# Add Control Bounds
################################################################################

function add_control_bound!(game_con::GameConstraintValues, u_max::AbstractVector, u_min::AbstractVector)
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	add_constraint!(game_con.control_conlist, ControlBoundConstraint(m,u_max=u_max,u_min=u_min), 1:N-1)
	con  = game_con.control_conlist.constraints[end]
	inds = game_con.control_conlist.inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.control_conval, conval)
	return nothing
end

################################################################################
# Add Circle Constraint
################################################################################

function add_circle_constraint!(game_con::GameConstraintValues, i::Int,
	xc::AbstractVector, yc::AbstractVector, radius::AbstractVector)
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	px = probsize.px

	add_constraint!(
		game_con.state_conlist[i],
		TrajectoryOptimization.CircleConstraint(n, xc, yc, radius, px[i][1], px[i][2]),
		2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end

function add_circle_constraint!(game_con::GameConstraintValues,
	xc::AbstractVector, yc::AbstractVector, radius::AbstractVector)
	p = game_con.probsize.p
	for i = 1:p
		add_circle_constraint!(game_con, i, xc, yc, radius)
	end
	return nothing
end

################################################################################
# Wall Constraint
################################################################################
abstract type AbstractWall end

mutable struct Wall <: AbstractWall
	p1::AbstractVector # initial point of the boundary
	p2::AbstractVector # final point of the boundary
	v::AbstractVector  # vector orthogonal to (p2 - p1) and indicating the forbiden halfspace
end

function add_wall_constraint!(game_con::GameConstraintValues, i::Int, walls::AbstractVector{Wall})
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	px = probsize.px

	n_wall = length(walls)
	T = eltype(walls[1].p1)
	x1 = SVector{n_wall,T}([wall.p1[1] for wall in walls])
	y1 = SVector{n_wall,T}([wall.p1[2] for wall in walls])
	x2 = SVector{n_wall,T}([wall.p2[1] for wall in walls])
	y2 = SVector{n_wall,T}([wall.p2[2] for wall in walls])
	xv = SVector{n_wall,T}([wall.v[1] for wall in walls])
	yv = SVector{n_wall,T}([wall.v[2] for wall in walls])

	add_constraint!(
		game_con.state_conlist[i],
		WallConstraint(n, x1, y1, x2, y2, xv, yv, px[i][1], px[i][2]),
		2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end

function add_wall_constraint!(game_con::GameConstraintValues, walls::AbstractVector{<:AbstractWall})
	p = game_con.probsize.p
	for i = 1:p
		add_wall_constraint!(game_con, i, walls)
	end
	return nothing
end

################################################################################
# Wall3D Constraint
################################################################################

mutable struct Wall3D <: AbstractWall
	p1::AbstractVector # initial point of the boundary
	p2::AbstractVector # middle point of the boundary
	p3::AbstractVector # final point of the boundary
	v::AbstractVector  # vector orthogonal to the plane defined by (p1,p2,p3) and indicating the forbiden halfspace
end

function add_wall_constraint!(game_con::GameConstraintValues, i::Int, walls::AbstractVector{Wall3D})
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	pz = probsize.pz

	n_wall = length(walls)
	T = eltype(walls[1].p1)
	x1 = SVector{n_wall,T}([wall.p1[1] for wall in walls])
	y1 = SVector{n_wall,T}([wall.p1[2] for wall in walls])
	z1 = SVector{n_wall,T}([wall.p1[3] for wall in walls])
	x2 = SVector{n_wall,T}([wall.p2[1] for wall in walls])
	y2 = SVector{n_wall,T}([wall.p2[2] for wall in walls])
	z2 = SVector{n_wall,T}([wall.p2[3] for wall in walls])
	x3 = SVector{n_wall,T}([wall.p3[1] for wall in walls])
	y3 = SVector{n_wall,T}([wall.p3[2] for wall in walls])
	z3 = SVector{n_wall,T}([wall.p3[3] for wall in walls])
	xv = SVector{n_wall,T}([wall.v[1] for wall in walls])
	yv = SVector{n_wall,T}([wall.v[2] for wall in walls])
	zv = SVector{n_wall,T}([wall.v[3] for wall in walls])

	add_constraint!(
		game_con.state_conlist[i],
		Wall3DConstraint(n,
			x1, y1, z1,
			x2, y2, z2,
			x3, y3, z3,
			xv, yv, zv, pz[i][1], pz[i][2], pz[i][3]),
		2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end
################################################################################
# CylinderWall Constraint
################################################################################

mutable struct CylinderWall <: AbstractWall
	p::AbstractVector # initial point of the boundary
	v::Symbol         # middle point of the boundary
	l::Real           # final point of the boundary
	r::Real           # vector orthogonal to the plane defined by (p1,p2,p3) and indicating the forbiden halfspace
end

function add_wall_constraint!(game_con::GameConstraintValues, i::Int, walls::AbstractVector{CylinderWall})
	probsize = game_con.probsize
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	pz = probsize.pz

	n_wall = length(walls)
	T = eltype(walls[1].p)
	p1 = SVector{n_wall,T}([wall.p[1] for wall in walls])
	p2 = SVector{n_wall,T}([wall.p[2] for wall in walls])
	p3 = SVector{n_wall,T}([wall.p[3] for wall in walls])
	v = SVector{n_wall,Symbol}([wall.v for wall in walls])
	l = SVector{n_wall,T}([wall.l for wall in walls])
	r = SVector{n_wall,T}([wall.r for wall in walls])

	add_constraint!(
		game_con.state_conlist[i],
		CylinderConstraint(n, p1, p2, p3, v, l, r, pz[i][1], pz[i][2], pz[i][3]),
		2:N)
	con  = game_con.state_conlist[i].constraints[end]
	inds = game_con.state_conlist[i].inds[end]
	conval = Altro.ALConVal(n,m,con,inds)
	push!(game_con.state_conval[i], conval)
	return nothing
end

################################################################################
# Helpers
################################################################################
function TrajectoryOptimization.cost_expansion!(conval::TrajectoryOptimization.AbstractConstraintValues)
    s = TrajectoryOptimization.sense(conval)
    for (i,k) in enumerate(conval.inds)
        TrajectoryOptimization.cost_expansion!(s, conval, i)
    end
    return nothing
end

function reset!(game_con::GameConstraintValues)
    reset_duals!(game_con)
    reset_penalties!(game_con)
    return nothing
end

function reset_duals!(game_con::GameConstraintValues)
    p = game_con.probsize.p
    for i = 1:p
        for conval in game_con.state_conval[i]
            Altro.reset_duals!(conval)
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        Altro.reset_duals!(conval)
    end
    return nothing
end

function reset_penalties!(game_con::GameConstraintValues)
    p = game_con.probsize.p
    for i = 1:p
        for conval in game_con.state_conval[i]
            Altro.reset_penalties!(conval)
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        Altro.reset_penalties!(conval)
    end
    return nothing
end

function penalty_update!(game_con::GameConstraintValues)
    p = game_con.probsize.p
    for i = 1:p
        for conval in game_con.state_conval[i]
            Altro.penalty_update!(conval)
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        Altro.penalty_update!(conval)
    end
    return nothing
end

function penalty_update!(con::Altro.ALConVal)
    for i in eachindex(con.μ)
        con.μ[i] .*= con.params.ϕ
    end
end

function dual_update!(game_con::GameConstraintValues)
    p = game_con.probsize.p
	α_dual = game_con.α_dual
	αx_dual = game_con.αx_dual
    for i = 1:p
        for conval in game_con.state_conval[i]
            # Altro.dual_update!(conval)
			dual_update!(conval, αx_dual[i])
        end
    end
    # Control constraints
    for conval in game_con.control_conval
		# Altro.dual_update!(conval)
		dual_update!(conval, α_dual)
    end
    return nothing
end

function evaluate!(game_con::GameConstraintValues, traj::Traj)
    p = game_con.probsize.p
    for i = 1:p
        for conval in game_con.state_conval[i]
            TrajectoryOptimization.evaluate!(conval, traj)
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        TrajectoryOptimization.evaluate!(conval, traj)
    end
    return nothing
end

import RobotDynamics.jacobian!
function jacobian!(game_con::GameConstraintValues, traj::Traj)
    p = game_con.probsize.p
    for i = 1:p
        for conval in game_con.state_conval[i]
            TrajectoryOptimization.jacobian!(conval, traj)
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        TrajectoryOptimization.jacobian!(conval, traj)
    end
    return nothing
end

function update_active_set!(game_con::GameConstraintValues, traj::Traj)
	evaluate!(game_con, traj)
	update_active_set!(game_con)
    return nothing
end


function update_active_set!(game_con::GameConstraintValues)
    p = game_con.probsize.p
    for i = 1:p
        for conval in game_con.state_conval[i]
            Altro.update_active_set!(conval, Altro.Val(game_con.active_set_tolerance))
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        Altro.update_active_set!(conval, Altro.Val(game_con.active_set_tolerance))
    end
    return nothing
end

################################################################################
# Dual Update
################################################################################

function dual_update!(conval::Altro.ALConVal, α_dual::T) where {T}
	c = conval.vals
	λ = conval.λ
	μ = conval.μ
	λ_max = conval.params.λ_max
	cone = TrajectoryOptimization.sense(conval.con)
	for i in eachindex(conval.inds)
		λ[i] .= dual_update(cone, SVector(λ[i]), SVector(c[i]), SVector(μ[i]), λ_max, α_dual)
	end
end

function dual_update(::Equality, λ, c, μ, λmax, α_dual)
	λbar = λ + α_dual*μ .* c
	return clamp.(λbar, -λmax, λmax)
end

function dual_update(::Inequality, λ, c, μ, λmax, α_dual)
 	λbar = λ + α_dual*μ .* c
	return clamp.(λbar, 0, λmax)  # project onto the dual cone via max(0,x)
end

function dual_update(cone::Altro.SecondOrderCone, λ, c, μ, λmax, α_dual)
	 λbar = λ - α_dual*μ .* c
	 return TrajectoryOptimization.projection(cone, λbar)  # project onto the dual cone
end
