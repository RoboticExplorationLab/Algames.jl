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
tf = 3.0  # final tim
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
dxf = @SVector [# p1   # p2   # p3
			   1.00,  1.00,  1.00, # x
			   0.00,  0.00,  0.00, # y
			   0.00,  0.00,  0.00, # θ
			   0.00,  0.00,  0.00, # v
			  ]

diag_Q = [SVector{n}([0.,  0.,  0.,
					  10., 0.,  0.,
					  1.,  0.,  0.,
					  1.,  0.,  0.]),
	      SVector{n}([0.,  0.,  0.,
		  			  0.,  10., 0.,
					  0.,  1.,  0.,
					  0.,  1.,  0.]),
		  SVector{n}([0.,  0.,  0.,
		  			  0.,  0.,  10.,
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
actors_types = [:car for i=1:p]
road_length = 6.0
road_width = 0.34
ramp_length = 3.2
ramp_angle = pi/12
ramp_merging_3_players_mpc_scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, actors_radii, actors_types)

# Progressive collision avodance radius
pr = 3 # 3 steps
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
actors_radii_prog = [[actor_radius + (j-1)*0.01 for i=1:p] for j=1:pr]

con_inds_prog = [j:j for j=1:pr]
con_inds_prog[end] = pr:N
# Create constraint sets
algames_conSet = ConstraintSet(n,m,N)
# Add collision avoidance constraints
for j = 1:pr
	add_collision_avoidance(algames_conSet, actors_radii_prog[j], px, p, con_inds_prog[j])
end
# Add scenario specific constraints (road boundaries)
con_inds = 2:N # Indices where the constraints will be applied
add_scenario_constraints(algames_conSet, ramp_merging_3_players_mpc_scenario, px, con_inds; constraint_type=:constraint)

# Create problems
algames_ramp_merging_3_players_mpc_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)

# Create solvers
algames_ramp_merging_3_players_mpc_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=0,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_ramp_merging_3_players_mpc_prob, algames_ramp_merging_3_players_mpc_opts)

# Create MPC solvers
state_noise = 5. * @SVector ([
	0.008, 0.008, 2*pi/72, 0.03, #+-50cm, +-50cm, +-25deg, +-12.5% per second
	0.008, 0.008, 2*pi/72, 0.03,
	0.008, 0.008, 2*pi/72, 0.03])
ramp_merging_3_players_mpc_opts = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=10,
	N_mpc=50,
	mpc_tf=100.0,
	min_δt=0.005,
	dxf=dxf,
	noise=state_noise)
algames_ramp_merging_3_players_mpc_solver = MPCGamesSolver(algames_solver, ramp_merging_3_players_mpc_opts)
# ilqgames_ramp_merging_3_players_mpc_solver = MPCGamesSolver(ilqgames_solver, dxf, opts_mpc)
# reset!(algames_ramp_merging_3_players_mpc_solver, reset_type=:full)
# reset!(ilqgames_ramp_merging_3_players_mpc_solver, reset_type=:full)
# solve!(algames_ramp_merging_3_players_mpc_solver; wait=false)
# resample!(algames_ramp_merging_3_players_mpc_solver)




#
#
# vis=AG.Visualizer()
# anim=AG.MeshCat.Animation()
# open(vis)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(mpc_solver, ramp_merging_3_players_mpc_scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=true)
#
# iter = mpc_solver.stats.iterations
# mean_solve_time = mean(mpc_solver.stats.solve_time[1:iter])
# update_freq = mpc_solver.stats.iterations/mpc_solver.stats.time
# std_solve_time = StatsBase.std(mpc_solver.stats.solve_time)
# largest_solve_time = maximum(mpc_solver.stats.solve_time)
#
#
# # To evaluate the update frequency of the MPC solver we average over many samples.
# samples = 100
# times = zeros(0)
# cmax = zeros(samples)
# mpc_solver.opts.log_level = AG.Logging.Warn
# for k = 1:samples
# 	@show k
# 	algames_solver = DirectGamesSolver(algames_prob, algames_opts)
# 	mpc_solver = MPCGamesSolver(algames_solver, dxf, opts_mpc)
#
#     reset!(mpc_solver, reset_type=:full)
#     solve!(mpc_solver, wait=false)
#     i = mpc_solver.stats.iterations
#     cmax[k] = maximum(mpc_solver.stats.cmax[1:i])
#     push!(times, mpc_solver.stats.solve_time[1:i]...)
# end
#
# # Average MPC frequency
# freq = length(times) / sum(times) # 174 Hz
# # Mean solve time
# mean_solve_time = sum(times) / length(times) #
# # Maximum constraint violation across samples.
# # This constraint violation should be comparable to the tolerance
# # we used in the solver.
# max_constraint_violation = maximum(cmax)
# solver_violation_tolerance = mpc_solver.solver.opts.constraint_tolerance
# # If the max_constraint_violation <= solver_violation_tolerance
# # this means that all the open-loop plans have converged to constraint
# # satisfaction during the MPC solve.
