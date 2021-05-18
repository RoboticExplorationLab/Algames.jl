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

function build_robot!(vis::Visualizer, model::DoubleIntegratorGame; name::Symbol=:Robot, α=0.8)
    p = model.p
	orange_mat, blue_mat, black_mat = get_material(;α=α)
    for i = 1:p
        setobject!(vis[name]["player$i"], Sphere(Point3f0(0.0), 0.08), orange_mat)
    end
    return nothing
end

function set_robot!(vis::Visualizer, model::DoubleIntegratorGame, s::AbstractVector; name::Symbol=:Robot)
    p = model.p
    pz = model.pz
    for i = 1:p
        x = s[pz[i][1:3]]
        settransform!(vis[name]["player$i"], MeshCat.Translation(x...))
    end
    return nothing
end
