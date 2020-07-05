using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES



# Define the dynamics model of the game.
struct TeamGame{N,M,P} <: AbstractGameModel
	n::Int   # Number of states
	n_::Int  # Number of states per player
	m::Int   # Number of controls
	m_::Int  # Number of controls per player
	p::Int   # Number of players per team
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px0::Vector{Int} # Indices of the passive player
	px::Vector{Vector{Int}} # Indices of the each player's states
end
TeamGame(;p::Int=2, d::Int=2) = TeamGame{d*(2*p+1),d*p,p}(
	# p = number of players per team
	# d = dimension of the integrator
	d*(2*p+1),
	2*d,
	d*p,
	d,
	p,
	[[(j-1)*2*p+i for j=1:d] for i=1:2*p],
	[1:2*d],
	[2*d .+ [(j-1)*2*p+i for j=1:2*d] for i=1:2*p],
	)
Base.size(model::DoubleIntegratorGame{N,M,P}) where {N,M,P} = N,M,
	model.pu,P # n,m,pu,p

@generated function TO.dynamics(model::TeamGame{N,M,P}, x, u) where {N,M,P}
	qd  = [:(x[$i]) for i=N-P*:N]
	qdd = [:(u[$i]) for i=1:M]
	return :(SVector{$N}($(), $(qd...), $(qdd...)))
end







# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Instantiate dynamics model
model =
model = DoubleIntegratorGame(p=2)
n,m,pu,p = size(model)
T = Float64

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [-0.50, -0,50,  0.00, -0.50, # x
			    0.10, -0,10,  0.00, -0.10, # y
                0.50,  0,40,  0.00,  0.40, # xdot
				0.00,  0,00,  0.00,  0.00] # ydot
xf = @SVector [ 0.50,  0.50,  0.00,  0.50, # x
			   -0.10,  0.10,  0.00,  0.10, # y
                0.40,  0.30,  0.00,  0.30, # xdot
				0.00,  0.10,  0.00,  0.10] # ydot

# Define a quadratic cost for each player
diag_Q = [SVector{n}([ 1.0, -0.1,  1.0, -0.1,  1.0, -0.1,  1.0, -0.1]), 	# Player 1 cost
	      SVector{n}([-0.1,  1.0, -0.1,  1.0, -0.1,  1.0, -0.1,  1.0])] 	# Player 2 cost
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N,checks=false) for i=1:p]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
ilqgames_conSet = ConstraintSet(n,m,N)
con_inds = 2:N

# # Add collision avoidance constraints
# add_collision_avoidance(algames_conSet, actors_radii, model.px, p, con_inds)
# add_collision_avoidance(ilqgames_conSet, actors_radii, model.px, p, con_inds)
# # # u_min = - SVector{m}(ones(m))
# # # u_max = + SVector{m}(ones(m))
# # # con = BoundConstraint(n,m,u_min=u_min,u_max=u_max)
# # # add_constraint!(algames_conSet, con, con_inds)
# # # add_constraint!(ilqgames_conSet, con, con_inds)

# Define the problem
algames_linear_quadratic_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)
ilqgames_linear_quadratic_prob = GameProblem(model, obj, xf, tf, constraints=ilqgames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    log_level=TO.Logging.Warn)
algames_linear_quadratic_solver = DirectGamesSolver(algames_linear_quadratic_prob, algames_opts)

ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn)
ilqgames_linear_quadratic_solver = PenaltyiLQGamesSolver(ilqgames_linear_quadratic_prob, ilqgames_opts)

solve!(algames_linear_quadratic_solver)
@time solve!(algames_linear_quadratic_solver)
solve!(ilqgames_linear_quadratic_solver)
@time solve!(ilqgames_linear_quadratic_solver)

visualize_trajectory_car(algames_linear_quadratic_solver)

# algames_linear_quadratic_solver.stats
# ilqgames_linear_quadratic_solver.stats

# X = TO.states(solver)
# U = TO.controls(solver)
# visualize_state(X)
# visualize_control(U,pu)
# visualize_trajectory_car(solver)
# visualize_collision_avoidance(solver)
# visualize_dynamics(solver)
# visualize_optimality_merit(solver)
# visualize_α(solver)
# visualize_cmax(solver)
