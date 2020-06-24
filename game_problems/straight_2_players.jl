using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Instantiate dynamics model
model = UnicycleGame(p=2)
n,m,pu,p = size(model)
T = Float64
px = model.px

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2
              -1.00, -1.00, # x
               0.10, -0.10, # y
			   0.00,  0.00, # θ
			   0.60,  0.60, # v
               ]
xf = @SVector [# p1   # p2
               1.00,  1.00, # x
               0.10, -0.05, # y
			   0.00,  0.00, # θ
			   0.60,  0.80, # v
              ]

diag_Q = [SVector{n}([0.,  0.,
					  10., 0.,
					  1.,  0.,
					  1.,  0.]),
	      SVector{n}([0.,  0.,
		  			  0.,  10.,
					  0.,  1.,
					  0.,  1.])]
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
actors_types = [:car for i=1:p]
road_length = 3.0
road_width = 0.37
obs_radius = 0.06
obs = [0.50, -0.17]
straight_2_players_scenario = StraightScenario(road_length, road_width,
	actors_radii, actors_types, obs_radius, obs)

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
ilqgames_conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
add_collision_avoidance(ilqgames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, straight_2_players_scenario, px, con_inds; constraint_type=:constraint)
add_scenario_constraints(ilqgames_conSet, straight_2_players_scenario, px, con_inds; constraint_type=:constraint)

algames_straight_2_players_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)
ilqgames_straight_2_players_prob = GameProblem(model, obj, xf, tf, constraints=ilqgames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Debug)
algames_straight_2_players_solver = DirectGamesSolver(algames_straight_2_players_prob, algames_opts)

ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn,
    )
ilqgames_straight_2_players_solver = PenaltyiLQGamesSolver(ilqgames_straight_2_players_prob, ilqgames_opts)
pen = ones(length(ilqgames_straight_2_players_solver.constraints))*1000.0
set_penalty!(ilqgames_straight_2_players_solver, pen);


# @time timing_solve(algames_straight_2_players_solver)
# @time timing_solve(ilqgames_straight_2_players_solver)
#
# visualize_trajectory_car(algames_straight_2_players_solver)
# visualize_trajectory_car(ilqgames_straight_2_players_solver)

#
# using Test
# using MeshCat
#
# prob = GameProblems.algames_straight_2_players_prob
#
# opts = DirectGamesSolverOptions(
#     iterations=10,
#     inner_iterations=20,
#     iterations_linesearch=10,
#     min_steps_per_iteration=1,
#     log_level=TO.Logging.Warn)
# solver = DirectGamesSolver(prob, opts)
# solve!(solver)

# vis = MeshCat.Visualizer()
# anim = MeshCat.Animation()
# # open(vis)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(solver, straight_2_players_scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=true)
#
# @test converged(solver)
