# Define the dynamics model of the game.
struct DoubleIntegratorGame{N,M,P} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px::Vector{Vector{Int}} # Indices of the each player's x and y positions
end
DoubleIntegratorGame(;p::Int=2, d::Int=2) = DoubleIntegratorGame{2*d*p,d*p,p}(
	# P = number of players
	# d = dimension of the integrator
	2*d*p,
	d*p,
	p,
	[[(j-1)*p+i for j=1:d] for i=1:p],
	[[(j-1)*p+i for j=1:d] for i=1:p])
Base.size(model::DoubleIntegratorGame{N,M,P}) where {N,M,P} = N,M,
	[[(j-1)*P+i for j=1:Int(M/P)] for i=1:P],P # n,m,pu,p

@generated function TO.dynamics(model::DoubleIntegratorGame{N,M,P}, x, u) where {N,M,P}
	qd  = [:(x[$i]) for i=M+1:N]
	qdd = [:(u[$i]) for i=1:M]
	return :(SVector{$N}($(qd...), $(qdd...)))
end
