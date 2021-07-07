# Define the dynamics model of the model.
struct DoubleIntegratorGame{N,M,P,SVu,SVx,SVz} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	ni::Vector{Int}  # Number of states for each player
	mi::Vector{Int}  # Number of controls for each player
	pu::SVu # Indices of the each player's controls
	px::SVx # Indices of the each player's x and y positions
	pz::SVz # Indices of the each player's states
end

function DoubleIntegratorGame(;p::Int=2, d::Int=2)
	# p = number of players
	# d = dimension of the integrator
	n = 2d*p
	m = d*p
	pu = [SVector{d,Int}([i + (j-1)*p for j=1:d]) for i=1:p]
	px = [SVector{2,Int}([i + (j-1)*p for j=1:2]) for i=1:p]
	pz = [SVector{2d,Int}([i + (j-1)*p for j=1:2d]) for i=1:p]
	TYPE = typeof.((pu,px,pz))
	ni = 2d*ones(Int,p)
	mi = d*ones(Int,p)
	return DoubleIntegratorGame{n,m,p,TYPE...}(n,m,p,ni,mi,pu,px,pz)
end

@generated function RobotDynamics.dynamics(model::DoubleIntegratorGame{N,M,P}, x, u) where {N,M,P}
	qd  = [:(x[$i]) for i=M+1:N]
	qdd = [:(u[$i]) for i=1:M]
	return :(SVector{$N}($(qd...), $(qdd...)))
end

dim(model::DoubleIntegratorGame) = model.mi[1]



# function build_robot!(vis::Visualizer, model::DoubleIntegratorGame; name::Symbol=:Robot, r_col=0.08, r_cost=0.5, α=1.0, r=0.15)
# 	r = convert(Float32, r)
#     p = model.p
#
# 	obj_path = joinpath(@__DIR__, "..", "mesh", "quadrotor", "drone.obj")
# 	mtl_path = joinpath(@__DIR__, "..", "mesh", "quadrotor", "drone.mtl")
# 	orange_mat, blue_mat, black_mat_col = get_material(;α=0.3)
# 	orange_mat, blue_mat, black_mat_cost = get_material(;α=0.1)
# 	obj = MeshFileObject(obj_path)
#     for i = 1:p
# 		setobject!(vis[name]["player$i"]["body"], obj)
# 		setobject!(vis[name]["player$i"]["collision"], GeometryBasics.Sphere(Point3f0(0.0), r_col), black_mat_col)
# 		setobject!(vis[name]["player$i"]["cost"], GeometryBasics.Sphere(Point3f0(0.0), r_cost), black_mat_cost)
# 		settransform!(vis[name]["player$i"]["body"], MeshCat.LinearMap(r*RotX(pi/2)))
#     end
#     return nothing
# end
#
# function set_robot!(vis::Visualizer, model::DoubleIntegratorGame, s::AbstractVector; name::Symbol=:Robot)
#     p = model.p
#     pz = model.pz
# 	d = dim(model)
#     for i = 1:p
# 		x = [s[pz[i][1:d]]; zeros(3-d)]
#         settransform!(vis[name]["player$i"], MeshCat.Translation(x...))
#     end
#     return nothing
# end
