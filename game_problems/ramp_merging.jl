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
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2   # p3
              -0.80, -1.00, -0.90, # x
              -0.05, -0.05, -0.30, # y
			   0.00,  0.00, pi/12, # θ
			   0.60,  0.60,  0.63, # v
               ]
xf = @SVector [# p1   # p2   # p3
               1.10,  0.70,  0.90, # x
              -0.05, -0.05, -0.05, # y
			   0.00,  0.00,  0.00, # θ
			   0.60,  0.60,  0.60, # v
              ]

diag_Q = [SVector{n}([0., 0., 0.,
					  1., 0., 0.,
					  1., 0., 0.,
					  1., 0., 0.]),
	      SVector{n}([0., 0., 0.,
		  			  0., 1., 0.,
					  0., 1., 0.,
					  0., 1., 0.]),
		  SVector{n}([0., 0., 0.,
		  			  0., 0., 1.,
					  0., 0., 1.,
					  0., 0., 1.])]
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
xs = SVector{n}(zeros(n))
us = SVector{m}(zeros(m))
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
actors_types = [:car, :car, :car]
road_length = 6.0
road_width = 0.30
ramp_length = 3.2
ramp_angle = pi/12
scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, actors_radii, actors_types)


# Create constraints
algames_conSet = ConstraintSet(n,m,N)
ilqgames_conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
add_collision_avoidance(ilqgames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, scenario, px, con_inds; constraint_type=:constraint)
add_scenario_constraints(ilqgames_conSet, scenario, px, con_inds; constraint_type=:constraint)

algames_ramp_merging = GameProblem(model, obj, algames_conSet, x0, xf, Z, N, tf)
ilqgames_ramp_merging = GameProblem(model, obj, ilqgames_conSet, x0, xf, Z, N, tf)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=0,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_ramp_merging, algames_opts)
ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn,
    )
ilqgames_solver = PenaltyiLQGamesSolver(ilqgames_ramp_merging, ilqgames_opts)
pen = ones(length(ilqgames_solver.constraints))*100.0
set_penalty!(ilqgames_solver, pen);

# @time timing_solve(algames_solver)
# @time timing_solve(ilqgames_solver)
#
# visualize_trajectory_car(algames_solver)
# visualize_trajectory_car(ilqgames_solver)









#
#
# model_iug = InertialUnicycleGame()
# x_iug = SVector{n}([i+0. for i=1:n])
# u_iug = SVector{m}([i+0. for i=1:m])
# x1_iug = TO.dynamics(model_iug, x_iug, u_iug)
#
# x1_iug_true = SVector{n}([
# 	4*cos(3),
# 	4*sin(3),
# 	1,
# 	2,
# 	8*cos(7),
# 	8*sin(7),
# 	3,
# 	4,
# 	12*cos(11),
# 	12*sin(11),
# 	5,
# 	6])
#
# x1_iug == x1_iug_true
#
# model_ug = UnicycleGame(p=3)
# x_ug = SVector{n}([i+0. for i=1:n])
# u_ug = SVector{m}([i+0. for i=1:m])
# x1_ug = TO.dynamics(model_ug, x_ug, u_ug)
#
#
# x1_ug_true = SVector{n}([
# 	10*cos(7),
# 	11*cos(8),
# 	12*cos(9),
# 	10*sin(7),
# 	11*sin(8),
# 	12*sin(9),
# 	1,
# 	2,
# 	3,
# 	4,
# 	5,
# 	6])
#
# x1_ug == x1_ug_true
# x1_ug == x1_iug
#
#
# size(model_ug)
