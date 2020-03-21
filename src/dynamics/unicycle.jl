# Define the dynamics model of the game.
struct UnicycleGame10{N,M,P} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
	p::Int  # Number of players
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px::Vector{Vector{Int}} # Indices of the each player's x and y positions
end
UnicycleGame10(;p::Int=2) = UnicycleGame10{4*p,2*p,p}(
	# P = number of players
	# d = dimension of the integrator
	4*p,
	2*p,
	p,
	[[(j-1)*2+i for j=1:2] for i=1:2],
	[[(j-1)*2+i for j=1:2] for i=1:2])
Base.size(::UnicycleGame10{N,M,P}) where {N,M,P} = N,M,
	[[(j-1)*Int(M/P)+i for j=1:Int(M/P)] for i=1:P],P # n,m,pu,p

@generated function TO.dynamics(model::UnicycleGame10{N,M,P}, x, u) where {N,M,P}
	xd  = [:(cos(x[$i])*x[$i+P]) for i=M+1:M+P]
	yd  = [:(sin(x[$i])*x[$i+P]) for i=M+1:M+P]
	qdd = [:(u[$i]) for i=1:M]
	return :(SVector{$N}($(xd...), $(yd...), $(qdd...)))
end
#
# # function TO.dynamics(model::UnicycleGame10, x, u) # Non memory allocating dynamics
# #     qd1 = @SVector [cos(x[3]), sin(x[3])]
# #     qd1 *= x[4]
# #     qd2 = @SVector [cos(x[7]), sin(x[7])]
# #     qd2 *= x[8]
# #     # qd3 = @SVector [cos(x[11]), sin(x[11])]
# #     # qd3 *= x[12]
# #     qdd1 = u[ @SVector [1,2] ]
# #     qdd2 = u[ @SVector [3,4] ]
# #     # qdd3 = u[ @SVector [5,6] ]
# #     return [qd1; qdd1; qd2; qdd2]#; qd3; qdd3]
# # end
#
#
# using BenchmarkTools
#
# using ALGAMES
# using LinearAlgebra
# using StaticArrays
# using TrajectoryOptimization
# const TO = TrajectoryOptimization
# const AG = ALGAMES
#
# # Discretization info
# tf = 3.0  # final time
# N = 41    # number of knot points
# dt = tf / (N-1) # time step duration
#
# # Instantiate dynamics model
model = UnicycleGame10(p=2)
# n,m,pu,p = size(model)
# T = Float64
#
# x = SVector{n}(zeros(n))
# u = SVector{m}(zeros(m))
#
# test_allocation(TO.dynamics, (model, x, u))
# test_allocation(TO.dynamics, (model, x, u))
# test_allocation(TO.dynamics, (model, x, u))
# test_allocation(TO.dynamics, (model, x, u))
# test_allocation(TO.dynamics, (model, x, u))
#
# @btime TO.dynamics($model, $x, $u)
#
