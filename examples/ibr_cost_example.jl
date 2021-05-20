################################################################################
# Iterative Best Response Example
################################################################################
# using Algames
using StaticArrays
using LinearAlgebra

T = Float64

# Define the dynamics of the system
p = 10 # Number of players
d = 2
model = DoubleIntegratorGame(p=p, d=d) # game with 3 players with unicycle dynamics
n = model.n
m = model.m

# Define the horizon of the problem
N = 20 # N time steps
dt = 0.1 # each step lasts 0.1 second
probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

# Define the objective of each player
# We use a LQR cost
Q = [Diagonal(10*ones(SVector{model.ni[i],T})) for i=1:p] # Quadratic state cost
R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
# Desrired state
xf = [SVector{model.ni[1],T}([ 2, 0.0,0,0]),
      SVector{model.ni[1],T}([ 2, 0.2,0,0]),
	  SVector{model.ni[1],T}([ 2, 0.4,0,0]),
	  SVector{model.ni[1],T}([ 2, 0.6,0,0]),
	  SVector{model.ni[1],T}([ 2, 0.8,0,0]),
	  SVector{model.ni[1],T}([ 2, 1.0,0,0]),
	  SVector{model.ni[1],T}([ 2, 1.2,0,0]),
	  SVector{model.ni[1],T}([ 2, 1.4,0,0]),
	  SVector{model.ni[1],T}([ 2, 1.6,0,0]),
	  SVector{model.ni[1],T}([ 2, 1.8,0,0]),
      ]
xf = xf[1:p]
# Desired control
uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
# Objectives of the game
game_obj = GameObjective(Q,R,xf,uf,N,model)
radius = 1.0*ones(p)
μ = 15.0*ones(p)
add_collision_cost!(game_obj, radius, μ)

# Define the constraints that each player must respect
game_con = GameConstraintValues(probsize)
# # Add collision avoidance
# radius = 0.08
# add_collision_avoidance!(game_con, radius)
# # Add control bounds
# u_max =  5*ones(SVector{m,T})
# u_min = -5*ones(SVector{m,T})
# add_control_bound!(game_con, u_max, u_min)
# # # Add state bounds for player 1
# x_max =  5*ones(SVector{n,T})
# x_min = -5*ones(SVector{n,T})
# add_state_bound!(game_con, 1, x_max, x_min)
# # Add wall constraint
# walls = [Wall([0.0,-0.4], [1.0,-0.4], [0.,-1.])]
# add_wall_constraint!(game_con, walls)
# # Add circle constraint
# xc = [1., 2., 3.]
# yc = [1., 2., 3.]
# radius = [0.1, 0.2, 0.3]
# add_circle_constraint!(game_con, xc, yc, radius)


# Define the initial state of the system
x0 = [
	[ 0.0,  0.0,  0.0,  0.0],
	[ 0.0,  0.2,  0.0,  0.0],
	[ 0.0,  0.4,  0.0,  0.0],
	[ 0.0,  0.6,  0.0,  0.0],
	[ 0.0,  0.8,  0.0,  0.0],
	[ 0.0,  1.0,  0.0,  0.0],
	[ 0.0,  1.2,  0.0,  0.0],
	[ 0.0,  1.4,  0.0,  0.0],
	[ 0.0,  1.6,  0.0,  0.0],
	[ 0.0,  1.8,  0.0,  0.0],
	]
x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))


# Define the Options of the solver
opts_alg = Options(Δ_min=1e-9, inner_print=false)
opts_ibr = Options(dual_reset=false, outer_iter=1, inner_iter=1, Δ_min=1e-5, inner_print=false)
# opts_ibr = Options(dual_reset=false, outer_iter=10, inner_iter=20, Δ_min=1e-5, inner_print=false)
ibr_opts = IBROptions(Δ_min=1e-9, ibr_iter=100, ordering=shuffle!([1,2,3,4,5,6,7,8,9,10]))
# ibr_opts = IBROptions(Δ_min=1e-9, ibr_iter=100, ordering=[1,2,3,4,5,6,7,8,9,10])
# ibr_opts = IBROptions(Δ_min=1e-9, ibr_iter=1, ordering=[1,2,3,4,5,6,7,8,9,10])

# Define the game problem
prob0_alg = GameProblem(N,dt,x0,model,opts_alg,game_obj,game_con)
prob0_ibr = GameProblem(N,dt,x0,model,opts_ibr,game_obj,game_con)
prob_alg = deepcopy(prob0_alg)
prob_ibr = deepcopy(prob0_ibr)

# cost.(prob_ibr.game_obj.obj)
# game_obj.obj

# Solve the problem
@elapsed newton_solve!(prob_alg)
@elapsed ibr_newton_solve!(prob_ibr, ibr_opts=ibr_opts)

residual!(prob_ibr, prob_ibr.pdtraj)
norm(prob_ibr.core.res, 1)/length(prob_ibr.core.res)
residual!(prob_alg, prob_alg.pdtraj)
norm(prob_alg.core.res, 1)/length(prob_alg.core.res)

# Visualize the Results
Algames.plot_traj!(prob_alg.model, prob_alg.pdtraj.pr)
Algames.plot_traj!(prob_ibr.model, prob_ibr.pdtraj.pr)
Algames.plot_violation!(prob_ibr.stats)
const Algames = Main
using Plots

prob_ibr.stats
prob_alg.stats
opts_ibr.dual_reset

for k = 1:N-1
	for i = 1:p
		prob_ibr.pdtraj.du[i][k] *= 0.0
	end
end

prob_ibr.pdtraj.du
#Lessons
	# - if there is no coupling at all (dynamics, cost, state constraints)
	# IBR marginally better, but should be the same since we exploit sparsity using a sparse linear solver.

	# - if there is cost coupling only
	# IBR converges to a different local NE
	# the solution depends on the players' ordering
	# the solution does not reflect the symmetry of information between players
	# some players having cost large advantages over

	# - if there is state constraints coupling
	#
