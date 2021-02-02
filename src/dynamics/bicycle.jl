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
