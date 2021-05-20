# Define the dynamics model of the model.
struct BicycleGame{N,M,P,LF,LR,SVu,SVx,SVz,T} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	ni::Vector{Int}  # Number of states for each player
	mi::Vector{Int}  # Number of controls for each player
	pu::SVu # Indices of the each player's controls
	px::SVx # Indices of the each player's x and y positions
	pz::SVz # Indices of the each player's states
	lf::T # Length from center of mass to front axle
	lr::T # Length from center of mass to rear axle
end

function BicycleGame(;p::Int=2, lf::T=0.05, lr::T=0.05) where {T}
	# p = number of players
	n = 4p
	m = 2p
	pu = [SVector{2,Int}([i + (j-1)*p for j=1:2]) for i=1:p]
	px = [SVector{2,Int}([i + (j-1)*p for j=1:2]) for i=1:p]
	pz = [SVector{4,Int}([i + (j-1)*p for j=1:4]) for i=1:p]
	TYPE = typeof.((pu,px,pz))
	ni = 4*ones(Int,p)
	mi = 2*ones(Int,p)
	return BicycleGame{n,m,p,lf,lr,TYPE...,T}(n,m,p,ni,mi,pu,px,pz,lf,lr)
end

@generated function RobotDynamics.dynamics(model::BicycleGame{N,M,P,lf,lr}, x, u) where {N,M,P,lf,lr}
	# https://archit-rstg.medium.com/two-to-four-bicycle-model-for-car-898063e87074
	# X = [x,y,v,ψ]
	# U = [a,δ]
	# β = atan(lr*tan(δ), lr+lf)
	# ̇X = [v*cos(β+ψ), v*cos(β+ψ), a, v*sin(β)/lr]

	L = :(lr+lf)
	xd  = [:(x[2P+$i]*cos(atan($lr*tan(u[P+$i]), $L) + x[3P+$i])) for i=1:P]
	yd  = [:(x[2P+$i]*sin(atan($lr*tan(u[P+$i]), $L) + x[3P+$i])) for i=1:P]
	vd  = [:(u[$i]) for i=1:P]
	ψd  = [:(x[2P+$i]*sin(atan($lr*tan(u[P+$i]), $L))/$lr) for i=1:P]
	return :(SVector{$N}($(xd...), $(yd...), $(vd...), $(ψd...)))
end

dim(model::BicycleGame) = 2


function build_robot!(vis::Visualizer, model::BicycleGame; name::Symbol=:Robot, r_col=0.08, r_cost=0.5, α=1.0, r=0.15)
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

function set_robot!(vis::Visualizer, model::BicycleGame, s::AbstractVector; name::Symbol=:Robot)
    p = model.p
    pz = model.pz
	d = dim(model)
    for i = 1:p
		x = [s[pz[i][1:d]]; zeros(3-d)]
		θ = s[pz[i][4]]
        settransform!(vis[name]["player$i"], MeshCat.compose(MeshCat.Translation(x...), MeshCat.LinearMap(RotZ(θ))))
    end
    return nothing
end
