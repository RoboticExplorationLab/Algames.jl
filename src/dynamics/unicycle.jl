# Define the dynamics model of the game.
struct UnicycleGame{N,M,P} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px::Vector{Vector{Int}} # Indices of the each player's x and y positions
end
UnicycleGame(;p::Int=2) = UnicycleGame{4*p,2*p,p}(
	# P = number of players
	# d = dimension of the integrator
	4*p,
	2*p,
	p,
	[[(j-1)*p+i for j=1:2] for i=1:p],
	[[(j-1)*p+i for j=1:2] for i=1:p])
Base.size(::UnicycleGame{N,M,P}) where {N,M,P} = N,M,
	[[(j-1)*P+i for j=1:2] for i=1:P],P # n,m,pu,p

@generated function TO.dynamics(model::UnicycleGame{N,M,P}, x, u) where {N,M,P}
	xd  = [:(cos(x[$i])*x[$i+P]) for i=M+1:M+P]
	yd  = [:(sin(x[$i])*x[$i+P]) for i=M+1:M+P]
	qdd = [:(u[$i]) for i=1:M]
	return :(SVector{$N}($(xd...), $(yd...), $(qdd...)))
end
