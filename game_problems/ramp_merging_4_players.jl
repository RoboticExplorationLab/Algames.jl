using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Instantiate dynamics model
model = UnicycleGame(p=4)
n,m,pu,p = size(model)
T = Float64
px = model.px

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2   # p3   # p4
              -0.80, -1.00, -0.90, -0.90,# x
              -0.10, -0.10, -0.35,  0.05,# y
			   0.00,  0.00, pi/12,  0.00,# θ
			   0.62,  0.60,  0.63,  0.79,# v
               ]
xf = @SVector [# p1   # p2   # p3   # p4
               1.10,  0.70,  0.90,  1.30,# x
              -0.10, -0.10, -0.10, -0.10,# y
			   0.00,  0.00,  0.00,  0.00,# θ
			   0.62,  0.60,  0.60,  0.79,# v
              ]

diag_Q = [SVector{n}([0.,  0.,  0.,  0.,
					  10., 0.,  0.,  0.,
					  1.,  0.,  0.,  0.,
					  1.,  0.,  0.,  0.]),
	      SVector{n}([0.,  0.,  0.,  0.,
		  			  0.,  10., 0.,  0.,
					  0.,  1.,  0.,  0.,
					  0.,  1.,  0.,  0.]),
		  SVector{n}([0.,  0.,  0.,  0.,
		  			  0.,  0.,  10., 0.,
					  0.,  0.,  1.,  0.,
					  0.,  0.,  1.,  0.]),
		  SVector{n}([0.,  0.,  0.,  0.,
					  0.,  0.,  0., 10.,
					  0.,  0.,  0.,  1.,
					  0.,  0.,  0.,  1.]),]
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
road_length = 6.0
road_width = 0.40
ramp_length = 3.2
ramp_angle = pi/12
ramp_merging_4_players_scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, actors_radii, actors_types)

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
ilqgames_conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
add_collision_avoidance(ilqgames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, ramp_merging_4_players_scenario, px, con_inds; constraint_type=:constraint)
add_scenario_constraints(ilqgames_conSet, ramp_merging_4_players_scenario, px, con_inds; constraint_type=:constraint)

algames_ramp_merging_4_players_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)
ilqgames_ramp_merging_4_players_prob = GameProblem(model, obj, xf, tf, constraints=ilqgames_conSet, x0=x0, N=N)

algames_ramp_merging_4_players_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=0,
    log_level=TO.Logging.Warn)
algames_ramp_merging_4_players_solver = DirectGamesSolver(algames_ramp_merging_4_players_prob, algames_ramp_merging_4_players_opts)

ilqgames_ramp_merging_4_players_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn,
    )
ilqgames_ramp_merging_4_players_solver = PenaltyiLQGamesSolver(ilqgames_ramp_merging_4_players_prob, ilqgames_ramp_merging_4_players_opts)
pen = ones(length(ilqgames_ramp_merging_4_players_solver.constraints))*100.0
set_penalty!(ilqgames_ramp_merging_4_players_solver, pen);

# @btime timing_solve(algames_ramp_merging_4_players_solver)
# @btime timing_solve(ilqgames_ramp_merging_4_players_solver)
#
# visualize_trajectory_car(algames_ramp_merging_4_players_solver)
# visualize_trajectory_car(ilqgames_ramp_merging_4_players_solver)
#
#
# vis, anim = animation(algames_ramp_merging_4_players_solver, scenario)
# vis, anim = animation(algames_ramp_merging_4_players_solver, scenario, vis=vis, anim=anim)
#
