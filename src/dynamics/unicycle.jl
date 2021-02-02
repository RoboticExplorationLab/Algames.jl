# Define the dynamics model of the model.
struct UnicycleGame{N,M,P,SVu,SVx,SVz} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	ni::Vector{Int}  # Number of states for each player
	mi::Vector{Int}  # Number of controls for each player
	pu::SVu # Indices of the each player's controls
	px::SVx # Indices of the each player's x and y positions
	pz::SVz # Indices of the each player's states
end

function UnicycleGame(;p::Int=2)
	# p = number of players
	n = 4p
	m = 2p
	pu = [SVector{2,Int}([i + (j-1)*p for j=1:2]) for i=1:p]
	px = [SVector{2,Int}([i + (j-1)*p for j=1:2]) for i=1:p]
	pz = [SVector{4,Int}([i + (j-1)*p for j=1:4]) for i=1:p]
	TYPE = typeof.((pu,px,pz))
	ni = 4*ones(Int,p)
	mi = 2*ones(Int,p)
	return UnicycleGame{n,m,p,TYPE...}(n,m,p,ni,mi,pu,px,pz)
end

@generated function RobotDynamics.dynamics(model::UnicycleGame{N,M,P}, x, u) where {N,M,P}
	xd  = [:(cos(x[$i])*x[$i+P]) for i=M+1:M+P]
	yd  = [:(sin(x[$i])*x[$i+P]) for i=M+1:M+P]
	qdd = [:(u[$i]) for i=1:M]
	return :(SVector{$N}($(xd...), $(yd...), $(qdd...)))
end
