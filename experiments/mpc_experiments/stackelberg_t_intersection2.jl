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
N = 61    # number of knot points
dt = tf / (N-1) # time step duration

diag_Q = [SVector{n}([0.,  0.,
					  1.,  0.,
					  1.,  0.,
					  10.,  0.]),
	      SVector{n}([0.,  0.,
		  			  0.,  1.,
					  0.,  1.,
					  0.,  10.])]
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
inflated_actors_radii = [2.6*actor_radius for i=1:p]
actors_types = [:car, :car]
top_road_length = 4.0
road_width = 0.44
bottom_road_length = 1.0
cross_width = 0.25
bound_radius = 0.05
lanes = [3,1]
t_intersection_2_players_scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length, cross_width, actors_radii,
		actors_types, bound_radius)

con_inds = 2:N # Indices where the constraints will be applied

algames_t_intersection_2_players_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
	μ_penalty=1.00,
    log_level=TO.Logging.Warn)

α_fading = 0.55

################################################################################
################################################################################
#####################           NASH EQUILIBRIUM          ######################
################################################################################
################################################################################



# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2
               0.11, -0.60, # x
              -0.38, -0.11, # y
			   pi/2,  0.00, # θ
			   0.20,  0.60, # v
               ]
xf = @SVector [# p1   # p2
              -1.30,  0.60, # x
               0.11, -0.11, # y
			     pi,  0.00, # θ
			   0.80,  0.60, # v
              ]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Create constraints
algames_conSet = ConstraintSet(n,m,N)

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, t_intersection_2_players_scenario,
	lanes, px, con_inds; constraint_type=:constraint)

algames_t_intersection_2_players_prob = GameProblem(
	model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)

algames_solver = DirectGamesSolver(
	algames_t_intersection_2_players_prob, algames_t_intersection_2_players_opts)

# add penalty constraints
add_collision_avoidance(algames_solver.penalty_constraints,
    inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)
reset!(algames_solver, reset_type=:full)

@time timing_solve(algames_solver)
visualize_trajectory_car(algames_solver)

using MeshCat
vis = MeshCat.Visualizer()
anim = MeshCat.Animation()
open(vis)
# sleep(1.0)
# Execute this line after the MeshCat tab is open
vis, anim = animation(algames_solver,
	t_intersection_2_players_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=false,
	camera_offset=false,
	camera_mvt=false,
	# α_fading=0.35
	)

################################################################################
################################################################################
#####################          STACKELBERG LEADER         ######################
################################################################################
################################################################################



# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2
		   0.11, -2.00, # x
		  -0.38, -0.11, # y
		   pi/2,  0.00, # θ
		   0.20,  0.00, # v
		   ]
xf = @SVector [# p1   # p2
		  -1.30,  0.60, # x
		   0.11, -0.11, # y
			 pi,  0.00, # θ
		   0.80,  0.00, # v
		  ]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Create constraints
algames_conSet = ConstraintSet(n,m,N)

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, t_intersection_2_players_scenario,
	lanes, px, con_inds; constraint_type=:constraint)

algames_t_intersection_2_players_prob = GameProblem(
	model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)

leader_solver = DirectGamesSolver(
	algames_t_intersection_2_players_prob, algames_t_intersection_2_players_opts)

# add penalty constraints
add_collision_avoidance(leader_solver.penalty_constraints,
	inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)
reset!(leader_solver, reset_type=:full)

@time timing_solve(leader_solver)
visualize_trajectory_car(leader_solver)

using MeshCat
vis_leader = MeshCat.Visualizer()
anim_leader = MeshCat.Animation()
open(vis_leader)
sleep(1.0)
# Execute this line after the MeshCat tab is open
vis_leader, anim_leader = animation(leader_solver,
	t_intersection_2_players_scenario;
	vis=vis_leader, anim=anim_leader,
	open_vis=false,
	display_actors=true,
	display_trajectory=false,
	camera_offset=false,
	camera_mvt=false,
	α_fading=0.35)



################################################################################
################################################################################
#####################         STACKELBERG FOLLOWER        ######################
################################################################################
################################################################################


# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2
		   0.11, -0.60, # x
		  -1.00, -0.11, # y
		   pi/2,  0.00, # θ
		   0.00,  0.60, # v
		   ]
xf = @SVector [# p1   # p2
		  -0.11,  0.60, # x
		  -1.00, -0.11, # y
		   pi/2,  0.00, # θ
		   0.00,  0.60, # v
		  ]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Create constraints
algames_conSet = ConstraintSet(n,m,N)

# Add the trajectory of the leader of the game as a constraint.
leader_states = deepcopy(TO.states(leader_solver))
leader_pos = [leader_states[k][[1,3]] for k = 1:N]
for k = 2:N
	leader_pos[k]
	radius = SVector{1}([actors_radii[1]+actors_radii[2]])
	xc =     SVector{1}([leader_pos[k][1]])
	yc =     SVector{1}([leader_pos[k][2]])
	col_con = CircleConstraint(n, xc, yc, radius, 2, 4)
	add_constraint!(algames_conSet, col_con, k:k)
end

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, t_intersection_2_players_scenario,
	lanes, px, con_inds; constraint_type=:constraint)

algames_t_intersection_2_players_prob = GameProblem(
	model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)

follower_solver = DirectGamesSolver(
	algames_t_intersection_2_players_prob, algames_t_intersection_2_players_opts)

# add penalty constraints
add_collision_avoidance(follower_solver.penalty_constraints,
	inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)
reset!(follower_solver, reset_type=:full)

@time timing_solve(follower_solver)
visualize_trajectory_car(follower_solver)

using MeshCat
vis_follower = MeshCat.Visualizer()
anim_follower = MeshCat.Animation()
open(vis_follower)
sleep(1.0)
# Execute this line after the MeshCat tab is open
vis_follower, anim_follower = animation(follower_solver,
	t_intersection_2_players_scenario;
	vis=vis_follower, anim=anim_follower,
	open_vis=false,
	display_actors=true,
	display_trajectory=false,
	camera_offset=false,
	camera_mvt=false,
	α_fading=0.30)
