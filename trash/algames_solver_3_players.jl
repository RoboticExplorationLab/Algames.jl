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
x0 = @SVector [# p1   # p2   # p3   # p4
	         -0.60, -1.20, -0.90, #-0.30, # x
		      0.08,  0.08, -0.31, #-0.08, # y
			  0.00,  0.00, pi/12, # 0.00, # θ
			  0.60,  0.60,  0.60, # 0.55, # v
]
xf = @SVector [# p1   # p2   # p3   # p4
			  1.10,  0.70,  1.30, # 1.15, # x
			 -0.08, -0.08,  0.08, # 0.08, # y
			  0.00,  0.00,  0.00, # 0.00, # θ
			  0.60,  0.60,  0.70, # 0.70, # v
              ]
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
					  0.,  0.,  1.]),]
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N,checks=false) for i=1:p]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
inflated_actors_radii = [1.2*actor_radius for i=1:p]
actors_types = [:car for i=1:p]
road_length = 34.20
road_width = 0.42
ramp_length = 17.2
ramp_angle = pi/12
ramp_merging_3_players_unicycle_penalty_scenario = MergingScenario(road_length,
	road_width, ramp_length, ramp_angle, actors_radii, actors_types)

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
con_inds = 1:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px,
	p, con_inds; constraint_type=:constraint)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, ramp_merging_3_players_unicycle_penalty_scenario,
	px, con_inds; constraint_type=:constraint)

# Add controls constraints
# u_lim = 0.12*ones(m)
u_lim = 0.50*ones(m)
control_bound = BoundConstraint(n,m,u_min=-u_lim,u_max=u_lim)
add_constraint!(algames_conSet, control_bound, 1:N-1)

algames_prob = GameProblem(model, obj, xf, tf,
	constraints=algames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
	optimality_constraint_tolerance=1e-2,
	μ_penalty=1.0,
    log_level=TO.Logging.Debug)
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

# add penalty constraints
add_collision_avoidance(algames_solver.penalty_constraints,
    inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)

reset!(algames_solver, reset_type=:full)
# algames_ramp_merging_2_players_unicycle_penalty_contraints = copy(algames_ramp_merging_2_players_unicycle_penalty_solver.penalty_constraints)

@time timing_solve(algames_solver)
visualize_trajectory_car(algames_solver)
visualize_control(algames_solver)

# using MeshCat
# vis = MeshCat.Visualizer()
# anim = MeshCat.Animation()
# open(vis)
# sleep(1.0)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(algames_solver,
# 	ramp_merging_3_players_unicycle_penalty_scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=true)
