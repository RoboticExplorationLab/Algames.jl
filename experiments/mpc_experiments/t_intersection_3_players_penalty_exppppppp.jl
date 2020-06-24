
using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Instantiate dynamics model
model = UnicycleGame(p=3)
n,m,pu,p = size(model)
T = Float64
px = model.px

# Discretization info
tf = 3.0  # final time
N = 21    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2   # p3
              -0.50,  1.40,  0.43, # x
              -0.15,  0.15, -0.30, # y
			   0.00,    pi,  pi/2, # x_dot
			   0.60,  0.60,  0.10, # y_dot
               ]
xf = @SVector [# p1   # p2   # p3
               1.30, -0.30,  0.43, # x
              -0.15,  0.15,  0.35, # y
			   0.00,    pi,  pi/2, # x_dot
			   0.60,  0.60,  0.30, # y_dot
              ]
# x0 = @SVector [
#                -0.50, -0.15,  0.00, 0.60, # player 1
#                 1.40,  0.15,  pi, 0.60,   # player 2
#                 0.43, -0.30,  pi/2, 0.10, # player 3
#                 ]
# xf = @SVector [
#                 1.30, -0.15,  0.00, 0.60, # player 1
#                -0.30,  0.15,  pi, 0.60,   # player 2
#                 0.43,  0.35,  pi/2, 0.30, # player 3
#                ]
diag_Q = [SVector{n}([0.,  0.,  0.,
					  1.,  0.,  0.,
					  1.,  0.,  0.,
					  1.,  0.,  0.]),
	      SVector{n}([0.,  0.,  0.,
		  			  0.,  1.,  0.,
					  0.,  1.,  0.,
					  0.,  1.,  0.]),
		  SVector{n}([0.,  0.,  0.,
		  			  0.,  0.,  1.,
					  0.,  0.,  1.,
					  0.,  0.,  1.])]
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
# inflated_actors_radii = [3.0*actor_radius for i=1:p] ####
inflated_actors_radii = [3.0*actor_radius for i=1:p]
actors_types = [:car, :car, :pedestrian]
top_road_length = 4.0
# road_width = 0.60 ###
road_width = 0.70
bottom_road_length = 1.0
cross_width = 0.25
bound_radius = 0.05
lanes = [1, 2, 5]
t_intersection_3_players_penalty_scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length, cross_width, actors_radii, actors_types, bound_radius)

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
con_inds = 1:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px,
	p, con_inds; constraint_type=:constraint)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, t_intersection_3_players_penalty_scenario,
	lanes, px, con_inds; constraint_type=:constraint)

algames_t_intersection_3_players_penalty_prob = GameProblem(model, obj, xf, tf,
	constraints=algames_conSet, x0=x0, N=N)

algames_t_intersection_3_players_penalty_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
	optimality_constraint_tolerance=1e-2,
	# μ_penalty=0.2,###
	μ_penalty=0.5,
    log_level=TO.Logging.Debug)
algames_t_intersection_3_players_penalty_solver =
	DirectGamesSolver(
	algames_t_intersection_3_players_penalty_prob,
	algames_t_intersection_3_players_penalty_opts)

# add penalty constraints
add_collision_avoidance(algames_t_intersection_3_players_penalty_solver.penalty_constraints,
    inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)

reset!(algames_t_intersection_3_players_penalty_solver, reset_type=:full)
algames_t_intersection_3_players_penalty_contraints = copy(algames_t_intersection_3_players_penalty_solver.penalty_constraints)

@time timing_solve(algames_t_intersection_3_players_penalty_solver)
visualize_trajectory_car(algames_t_intersection_3_players_penalty_solver)


using MeshCat
vis = MeshCat.Visualizer()
anim = MeshCat.Animation()
open(vis)
sleep(1.0)
# Execute this line after the MeshCat tab is open
vis, anim = animation(algames_t_intersection_3_players_penalty_solver,
	t_intersection_3_players_penalty_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=false,
	camera_offset=false,
	α_fading=1.00)
