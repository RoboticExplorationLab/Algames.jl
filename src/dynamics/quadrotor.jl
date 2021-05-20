
# Define the dynamics model of the model.
struct QuadrotorGame{N,M,P,SVu,SVx,SVz,T} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	ni::Vector{Int}  # Number of states for each player
	mi::Vector{Int}  # Number of controls for each player
	pu::SVu # Indices of the each player's controls
	px::SVx # Indices of the each player's x and y positions
	pz::SVz # Indices of the each player's states
	mass::T # mass of the quadrotor, in kg (default = 0.5)
    J::Diagonal{T,SVector{3,T}} # inertia of the quadrotor, in kg⋅m² (default = `Diagonal([0.0023, 0.0023, 0.004])`)
    Jinv::Diagonal{T,SVector{3,T}} #
    gravity::SVector{3,T} # gravity vector, in kg/m² (default = [0,0,-9.81])
    motor_dist::T # distance between the motors, in m (default = 0.1750)
    kf::T # motor torque constant (default = 0.0245)
    km::T # motor force constant (default = 1.0)
end

function QuadrotorGame(;p::Int=2, mass::T=0.5) where {T}
	@assert p <= 4 # if p > 4, you need to define the dynamics function associated with it.
	J = Diagonal(SVector{3,T}([0.0023, 0.0023, 0.004]))
	J_inv = Diagonal(SVector{3,T}(1 ./ [0.0023, 0.0023, 0.004]))
	gravity = SVector{3,T}([0,0,-9.81])
	motor_dist = 0.1750
	# kf = 0.0245
	# kf = 0.245
	kf = 1.245
	km = 1.0

	# p = number of players
	n = 12p
	m = 4p
	pu = [SVector{4,Int}([i + (j-1)*p for j=1:4]) for i=1:p]
	px = [SVector{2,Int}([i + (j-1)*p for j=1:2]) for i=1:p]
	pz = [SVector{12,Int}([i + (j-1)*p for j=1:12]) for i=1:p]
	TYPE = typeof.((pu,px,pz))
	ni = 12*ones(Int,p)
	mi = 4*ones(Int,p)
	return QuadrotorGame{n,m,p,TYPE...,T}(
		n,m,p,ni,mi,pu,px,pz,
		mass, J, J_inv, gravity, motor_dist, kf, km)
end


inertia(model::QuadrotorGame) = model.J
inertia_inv(model::QuadrotorGame) = model.Jinv
mass(model::QuadrotorGame) = model.mass

function forces(model::QuadrotorGame{N,M,P,SVu,SVx,SVz,T}, x, u, i::Int) where {N,M,P,SVu,SVx,SVz,T}
	R = MRP(x[3P+i], x[4P+i], x[5P+i]) # MRP indices for player i
    kf = model.kf
    g = model.gravity
    m = model.mass

    w1 = u[0P+i]
    w2 = u[1P+i]
    w3 = u[2P+i]
    w4 = u[3P+i]

    F1 = max(0,kf*w1);
    F2 = max(0,kf*w2);
    F3 = max(0,kf*w3);
    F4 = max(0,kf*w4);
    F = @SVector [0., 0., F1+F2+F3+F4] #total rotor force in body frame

    f = m*g + R*F # forces in world frame
    return f
end

function moments(model::QuadrotorGame{N,M,P,SVu,SVx,SVz,T}, x, u, i::Int) where {N,M,P,SVu,SVx,SVz,T}
    kf, km = model.kf, model.km
    L = model.motor_dist

    w1 = u[0P+i]
    w2 = u[1P+i]
    w3 = u[2P+i]
    w4 = u[3P+i]

    F1 = max(0,kf*w1);
    F2 = max(0,kf*w2);
    F3 = max(0,kf*w3);
    F4 = max(0,kf*w4);

    M1 = km*w1;
    M2 = km*w2;
    M3 = km*w3;
    M4 = km*w4;
    tau = @SVector [L*(F2-F4), L*(F3-F1), (M1-M2+M3-M4)] #total rotor torque in body frame
    return tau
end

function wrenches(model::QuadrotorGame, x::SVector, u::SVector, i::Int)
    F = forces(model, x, u, i)
    M = moments(model, x, u, i)
    return F, M
end

function RobotDynamics.dynamics(model::QuadrotorGame{N,M,P,SVu,SVx,SVz,T},
		x::SVector, u::SVector, i::Int) where {N,M,P,SVu,SVx,SVz,T}
	#
	q = @SVector [x[3P+i], x[4P+i], x[5P+i]]
	v = @SVector [x[6P+i], x[7P+i], x[8P+i]]
	ω = @SVector [x[9P+i], x[10P+i], x[11P+i]]
	R = MRP(q)

    # Original dynamics
    F,τ = wrenches(model, x, u, i)
    m = mass(model)
    J = inertia(model)
    Jinv = inertia_inv(model)

    xdot = v
    qdot = Rotations.kinematics(R,ω)
    vdot = F ./ m
    ωdot = Jinv*(τ - ω × (J*ω))

	return xdot, qdot, vdot, ωdot
end

function RobotDynamics.dynamics(model::QuadrotorGame{N,M,1,SVu,SVx,SVz,T},
		x::SVector, u::SVector) where {N,M,SVu,SVx,SVz,T}

	xdot1, qdot1, vdot1, ωdot1 = RobotDynamics.dynamics(model, x, u, 1)
	return @SVector [
		xdot1[1],
		xdot1[2],
		xdot1[3],
		qdot1[1],
		qdot1[2],
		qdot1[3],
		vdot1[1],
		vdot1[2],
		vdot1[3],
		ωdot1[1],
		ωdot1[2],
		ωdot1[3],
		]
end

function RobotDynamics.dynamics(model::QuadrotorGame{N,M,2,SVu,SVx,SVz,T},
		x::SVector, u::SVector) where {N,M,SVu,SVx,SVz,T}

	xdot1, qdot1, vdot1, ωdot1 = RobotDynamics.dynamics(model, x, u, 1)
	xdot2, qdot2, vdot2, ωdot2 = RobotDynamics.dynamics(model, x, u, 2)
	return @SVector [
		xdot1[1], xdot2[1],
		xdot1[2], xdot2[2],
		xdot1[3], xdot2[3],
		qdot1[1], qdot2[1],
		qdot1[2], qdot2[2],
		qdot1[3], qdot2[3],
		vdot1[1], vdot2[1],
		vdot1[2], vdot2[2],
		vdot1[3], vdot2[3],
		ωdot1[1], ωdot2[1],
		ωdot1[2], ωdot2[2],
		ωdot1[3], ωdot2[3],
		]
end

function RobotDynamics.dynamics(model::QuadrotorGame{N,M,3,SVu,SVx,SVz,T},
		x::SVector, u::SVector) where {N,M,SVu,SVx,SVz,T}

	xdot1, qdot1, vdot1, ωdot1 = RobotDynamics.dynamics(model, x, u, 1)
	xdot2, qdot2, vdot2, ωdot2 = RobotDynamics.dynamics(model, x, u, 2)
	xdot3, qdot3, vdot3, ωdot3 = RobotDynamics.dynamics(model, x, u, 3)
	return @SVector [
		xdot1[1], xdot2[1], xdot3[1],
		xdot1[2], xdot2[2], xdot3[2],
		xdot1[3], xdot2[3], xdot3[3],
		qdot1[1], qdot2[1], qdot3[1],
		qdot1[2], qdot2[2], qdot3[2],
		qdot1[3], qdot2[3], qdot3[3],
		vdot1[1], vdot2[1], vdot3[1],
		vdot1[2], vdot2[2], vdot3[2],
		vdot1[3], vdot2[3], vdot3[3],
		ωdot1[1], ωdot2[1], ωdot3[1],
		ωdot1[2], ωdot2[2], ωdot3[2],
		ωdot1[3], ωdot2[3], ωdot3[3],
		]
end

function RobotDynamics.dynamics(model::QuadrotorGame{N,M,4,SVu,SVx,SVz,T},
		x::SVector, u::SVector) where {N,M,SVu,SVx,SVz,T}

	xdot1, qdot1, vdot1, ωdot1 = RobotDynamics.dynamics(model, x, u, 1)
	xdot2, qdot2, vdot2, ωdot2 = RobotDynamics.dynamics(model, x, u, 2)
	xdot3, qdot3, vdot3, ωdot3 = RobotDynamics.dynamics(model, x, u, 3)
	xdot4, qdot4, vdot4, ωdot4 = RobotDynamics.dynamics(model, x, u, 4)
	return @SVector [
		xdot1[1], xdot2[1], xdot3[1], xdot4[1],
		xdot1[2], xdot2[2], xdot3[2], xdot4[2],
		xdot1[3], xdot2[3], xdot3[3], xdot4[3],
		qdot1[1], qdot2[1], qdot3[1], qdot4[1],
		qdot1[2], qdot2[2], qdot3[2], qdot4[2],
		qdot1[3], qdot2[3], qdot3[3], qdot4[3],
		vdot1[1], vdot2[1], vdot3[1], vdot4[1],
		vdot1[2], vdot2[2], vdot3[2], vdot4[2],
		vdot1[3], vdot2[3], vdot3[3], vdot4[3],
		ωdot1[1], ωdot2[1], ωdot3[1], ωdot4[1],
		ωdot1[2], ωdot2[2], ωdot3[2], ωdot4[2],
		ωdot1[3], ωdot2[3], ωdot3[3], ωdot4[3],
		]
end

dim(model::QuadrotorGame) = 3



function build_robot!(vis::Visualizer, model::QuadrotorGame; name::Symbol=:Robot, r_col=0.08, r_cost=0.5, α=1.0, r=0.15)
	r = convert(Float32, r)
    p = model.p

	obj_path = joinpath(@__DIR__, "..", "mesh", "quadrotor", "drone.obj")
	mtl_path = joinpath(@__DIR__, "..", "mesh", "quadrotor", "drone.mtl")
	orange_mat, blue_mat, black_mat_col = get_material(;α=0.3)
	orange_mat, blue_mat, black_mat_cost = get_material(;α=0.1)
	obj = MeshFileObject(obj_path)
    for i = 1:p
		setobject!(vis[name]["player$i"]["body"], obj)
		setobject!(vis[name]["player$i"]["collision"], GeometryBasics.Sphere(Point3f0(0.0), r_col), black_mat_col)
		setobject!(vis[name]["player$i"]["cost"], GeometryBasics.Sphere(Point3f0(0.0), r_cost), black_mat_cost)
		settransform!(vis[name]["player$i"]["body"], MeshCat.LinearMap(r*RotX(pi/2)))
    end
    return nothing
end

function set_robot!(vis::Visualizer, model::QuadrotorGame, s::AbstractVector; name::Symbol=:Robot)
    p = model.p
    pz = model.pz
    for i = 1:p
		x = s[pz[i][1:3]]
		R = MRP(s[pz[i][4:6]]...)
        settransform!(vis[name]["player$i"], MeshCat.compose(MeshCat.Translation(x...), MeshCat.LinearMap(R)))
    end
    return nothing
end
