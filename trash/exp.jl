using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Instantiate dynamics model
model = DoubleIntegratorGame(p=2)
n,m,pu,p = size(model)
T = Float64

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [-0.50, -0.50, # x
			    0.10, -0.10, # y
                0.50,  0.40, # xdot
				0.00,  0.00] # ydot
xf = @SVector [ 0.50,  0.50, # x
			   -0.10,  0.10, # y
                0.40,  0.30, # xdot
				0.00,  0.10] # ydot

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
algames_lq_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)
ilqgames_lq_prob = GameProblem(model, obj, xf, tf, constraints=ilqgames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    log_level=TO.Logging.Debug)
algames_lq_solver = DirectGamesSolver(algames_lq_prob, algames_opts)

ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Debug)
ilqgames_lq_solver = PenaltyiLQGamesSolver(ilqgames_lq_prob, ilqgames_opts)

solve!(algames_lq_solver)
@time solve!(algames_lq_solver)
solve!(ilqgames_lq_solver)
@time solve!(ilqgames_lq_solver)

# algames_lq_solver.stats
# ilqgames_lq_solver.stats

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




#
# for i in eachindex(algames_lq_solver.constraints.constraints)
# 	con = algames_lq_solver.constraints.constraints[i]
# 	@show con
# end
#
# for con in algames_lq_solver.constraints.constraints
# 	@show con
# end
#
# for con in algames_lq_solver.constraints
# 	@show con
# end
#
# if !isempty(algames_lq_solver.constraints.constraints)
# 	cmax = max(maximum(algames_lq_solver.constraints.c_max),
# 		maximum(algames_lq_solver.dyn_constraints.c_max))
# end
#
# #######
#
# for i in eachindex(ilqgames_lq_solver.constraints.constraints)
# 	con = ilqgames_lq_solver.constraints.constraints[i]
# 	@show con
# end
#
# for con in ilqgames_lq_solver.constraints.constraints
# 	@show con
# end
#
# for con in ilqgames_lq_solver.constraints
# 	@show con
# end
#
# if !isempty(ilqgames_lq_solver.constraints.constraints)
# 	@show 1
# end




algames_mpc_solver = GameProblems.algames_ramp_merging_3_players_mpc_solver
reset!(algames_mpc_solver, reset_type=:full)
algames_mpc_solver.solver.opts.constraint_tolerance = 1e-3
algames_mpc_solver.solver.opts.optimality_constraint_tolerance = 1e-2
algames_mpc_solver.solver.opts.log_level = Logging.Debug
solve!(algames_mpc_solver; wait=false)
resample!(algames_mpc_solver)




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
scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, actors_radii, actors_types)

# Progressive collision avodance radius
pr = 3 # 3 steps
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
actors_radii_prog = [[actor_radius + (j-1)*0.01 for i=1:p] for j=1:pr]

con_inds_prog = [j:j for j=1:pr]
con_inds_prog[end] = pr:N
# Create constraint sets
algames_conSet = ConstraintSet(n,m,N)
# ilqgames_conSet = ConstraintSet(n,m,N)
# Add collision avoidance constraints
for j = 1:pr
	add_collision_avoidance(algames_conSet, actors_radii_prog[j], px, p, con_inds_prog[j])
	# add_collision_avoidance(ilqgames_conSet, actors_radii_prog[j], px, p, con_inds_prog[j])
end
# Add scenario specific constraints (road boundaries)
con_inds = 2:N # Indices where the constraints will be applied
add_scenario_constraints(algames_conSet, scenario, px, con_inds; constraint_type=:constraint)
# add_scenario_constraints(ilqgames_conSet, scenario, px, con_inds; constraint_type=:constraint)

# Create problems
algames_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)
# ilqgames_prob = GameProblem(model, obj, xf, tf, constraints=ilqgames_conSet, x0=x0, N=N)

# Create solvers
algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=0,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

# ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
#     iterations=200,
#     gradient_norm_tolerance=1e-2,
#     cost_tolerance=1e-4,
#     line_search_lower_bound=0.0,
#     line_search_upper_bound=0.05,
#     log_level=TO.Logging.Warn,
#     )
# ilqgames_solver = PenaltyiLQGamesSolver(ilqgames_prob, ilqgames_opts)
# pen = ones(length(ilqgames_solver.constraints))*1000.0
# set_penalty!(ilqgames_solver, pen);

# Create MPC solvers
state_noise = 5. * @SVector ([
	0.008, 0.008, 2*pi/72, 0.03, #+-50cm, +-50cm, +-25deg, +-12.5% per second
	0.008, 0.008, 2*pi/72, 0.03,
	0.008, 0.008, 2*pi/72, 0.03])
opts_mpc = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=10,
	N_mpc=50,
	mpc_tf=100.0,
	min_δt=0.005,
	noise=state_noise)
algames_ramp_merging_3_players_mpc_solver = MPCGamesSolver(algames_solver, dxf, opts_mpc)
# ilqgames_ramp_merging_3_players_mpc_solver = MPCGamesSolver(ilqgames_solver, dxf, opts_mpc)
# reset!(algames_ramp_merging_3_players_mpc_solver, reset_type=:full)
# reset!(ilqgames_ramp_merging_3_players_mpc_solver, reset_type=:full)
# solve!(algames_ramp_merging_3_players_mpc_solver; wait=false)
# resample!(algames_ramp_merging_3_players_mpc_solver)



















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
              -0.80, -0.90, # x
              -0.05, -0.30, # y
			   0.00, pi/12, # θ
			   0.60,  0.63, # v
               ]
xf = @SVector [# p1   # p2
               1.10,  0.90, # x
              -0.05, -0.05, # y
			   0.00,  0.00, # θ
			   0.60,  0.60, # v
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
road_length = 6.0
road_width = 0.34
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

algames_ramp_merging_2_players_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)
ilqgames_ramp_merging_2_players_prob = GameProblem(model, obj, xf, tf, constraints=ilqgames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=40,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Debug)
algames_ramp_merging_2_players_solver = DirectGamesSolver(algames_ramp_merging_2_players_prob, algames_opts)

ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn,
    )
ilqgames_ramp_merging_2_players_solver = PenaltyiLQGamesSolver(ilqgames_ramp_merging_2_players_prob, ilqgames_opts)
pen = ones(length(ilqgames_ramp_merging_2_players_solver.constraints))*1000.0
set_penalty!(ilqgames_ramp_merging_2_players_solver, pen);

for con in algames_ramp_merging_2_players_solver.constraints.constraints
	con.params.μ_max = 1e40
end


@time timing_solve(algames_ramp_merging_2_players_solver)
# @time timing_solve(ilqgames_ramp_merging_2_players_solver)
#
# visualize_trajectory_car(algames_ramp_merging_2_players_solver)
# visualize_trajectory_car(ilqgames_ramp_merging_2_players_solver)
# converged(algames_ramp_merging_2_players_solver)
# X = TO.states(algames_ramp_merging_2_players_solver)
# U = TO.controls(algames_ramp_merging_2_players_solver)
# visualize_state(X)
# visualize_control(U,pu)
# visualize_trajectory_car(algames_ramp_merging_2_players_solver)
# visualize_collision_avoidance(algames_ramp_merging_2_players_solver)
# visualize_dynamics(algames_ramp_merging_2_players_solver)
# visualize_optimality_merit(algames_ramp_merging_2_players_solver)
# visualize_α(algames_ramp_merging_2_players_solver)
# visualize_cmax(algames_ramp_merging_2_players_solver)
#
# con = algames_ramp_merging_2_players_solver.constraints.constraints[1]
# con.μ
# solver = algames_ramp_merging_2_players_solver
# cost_expansion(solver.C, solver.obj, solver.Z, solver.model.pu, solver.model.p)
# TO.evaluate!(solver.constraints, solver.Z)
# TO.jacobian!(solver.constraints, solver.Z)
# TO.update_active_set!(solver.constraints,solver.Z)
# TO.evaluate!(solver.dyn_constraints, solver.Z)
# TO.discrete_jacobian!(solver.∇F, solver.model, solver.Z)
#
# update_g_!(solver)
# using StatsBase
# maximum(abs.(solver.g_))
#








using Test

atol = 1e-14
# Solving a linear quadratic game with ALGAMES
#  - Should be done in one outer iteration
#  - Should be done in 2 inner iterations (one to take a step the second one to check we are at the optimal point)
#  - Should converge perfectly: optimalty constraint ≈ 0.0
algames_prob = GameProblems.algames_linear_quadratic_prob
algames_opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

reset!(algames_solver, reset_type=:full)
solve!(algames_solver)
iter = algames_solver.stats.iterations
inner_iter = algames_solver.stats.iterations_inner[iter]
@test iter == 1
@test inner_iter == 2
@test isapprox(0., algames_solver.stats.optimality_merit[iter][inner_iter], atol=atol)


# Solving a linear quadratic game with ilQGames
#  - Should be done in 2 inner iterations (one to take a step the second one to check we are at the optimal point)
#  - Should converge perfectly: optimalty constraint ≈ 0.0
ilqgames_prob = GameProblems.ilqgames_linear_quadratic_prob
ilqgames_opts = PenaltyiLQGamesSolverOptions(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn)
ilqgames_solver = PenaltyiLQGamesSolver(ilqgames_prob, ilqgames_opts)

reset!(ilqgames_solver, reset_type=:full)
ilqgames_solver.opts.line_search_upper_bound=Inf
solve!(ilqgames_solver)
iter = ilqgames_solver.stats.iterations
@test iter == 2
@test isapprox(0., ilqgames_solver.stats.gradient[iter], atol=atol)














algames_prob = GameProblems.algames_ramp_merging_3_players_prob
opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_prob, opts)

ilqgames_prob = GameProblems.ilqgames_ramp_merging_3_players_prob
opts = PenaltyiLQGamesSolverOptions(
    iterations=600,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.02,
    log_level=TO.Logging.Warn)
pen = ones(length(ilqgames_prob.constraints))*1000.0
ilqgames_solver = PenaltyiLQGamesSolver(ilqgames_prob, opts)
set_penalty!(ilqgames_solver, pen)

num_samples = 10
state_noise = @SVector [ # Uniform noise around x0
    0.06, 0.06, 2*pi/72, 0.05,
    0.06, 0.06, 2*pi/72, 0.05,
    0.06, 0.06, 2*pi/72, 0.05]

opts_monte_carlo = MonteCarloSamplerOptions(
    noise=state_noise, # noise added to the initial state
    iterations=num_samples) # number of Monte Carlo samples

algames_sampler = MonteCarloSampler(algames_solver, opts_monte_carlo)
ilqgames_sampler = MonteCarloSampler(ilqgames_solver, opts_monte_carlo)

monte_carlo_sampling(algames_sampler)
monte_carlo_sampling(ilqgames_sampler)

@test num_converged(algames_sampler) == num_samples
@test num_converged(ilqgames_sampler) == num_samples








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
scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, actors_radii, actors_types)

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
add_scenario_constraints(algames_conSet, scenario, px, con_inds; constraint_type=:constraint)

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
# vis, anim = animation(mpc_solver, scenario;
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


using Test

algames_prob = GameProblems.algames_ramp_merging_3_players_mpc_prob
algames_opts = GameProblems.algames_ramp_merging_3_players_mpc_opts
mpc_opts = GameProblems.ramp_merging_3_players_mpc_opts

# Create solver
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

# We solve the problem once before to avoid a substantial part of
# Julia's precompilation time for the MPC test.
solve!(algames_solver)
reset!(algames_solver, reset_type=:full)

# Create MPC solver
mpc_solver = MPCGamesSolver(algames_solver, mpc_opts)
reset!(mpc_solver, reset_type=:full)
# mpc_solver.opts.log_level = Logging.Warn
mpc_solver.solver.opts.log_level = Logging.Warn
solve!(mpc_solver; wait=true)
resample!(mpc_solver)
# We should be able to solve the 10 MPC steps in less than mpc_solver.opts.mpc_tf = 100s.
@test mpc_solver.stats.time <= mpc_solver.opts.mpc_tf



























prob = GameProblems.algames_straight_2_players_prob

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)
solver = DirectGamesSolver(prob, opts)
solve!(solver)

@test converged(solver)



using Test

prob = GameProblems.algames_straight_2_players_prob

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)
solver = DirectGamesSolver(prob, opts)
solve!(solver)

visualize_state(solver)
visualize_control(solver)
visualize_trajectory_car(solver)
visualize_collision_avoidance(solver)
visualize_circle_collision(solver)
visualize_boundary_collision(solver)
visualize_dynamics(solver)
visualize_optimality_merit(solver)
visualize_H_cond(solver)
visualize_α(solver)
visualize_cmax(solver)

@test converged(solver)








using Test
using MeshCat

prob = GameProblems.algames_straight_2_players_prob
scenario = GameProblems.straight_2_players_scenario

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)
solver = DirectGamesSolver(prob, opts)
solve!(solver)

vis = MeshCat.Visualizer()
anim = MeshCat.Animation()
# open(vis)
# Execute this line after the MeshCat tab is open
vis, anim = animation(solver, scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

@test converged(solver)






algames_prob = GameProblems.algames_ramp_merging_3_players_mpc_prob
algames_opts = GameProblems.algames_ramp_merging_3_players_mpc_opts
mpc_opts = GameProblems.ramp_merging_3_players_mpc_opts

# Create solver
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

# We solve the problem once before to avoid a substantial part of
# Julia's precompilation time for the MPC test.
solve!(algames_solver)
reset!(algames_solver, reset_type=:full)

# Create MPC solver.
mpc_solver = MPCGamesSolver(algames_solver, mpc_opts)
mpc_solver.solver.opts.log_level = Logging.Warn
mpc_solver.opts.max_δt = 0.05
reset!(mpc_solver, reset_type=:full)
mpc_solver.opts.iterations = 10
solve!(mpc_solver; wait=true)
resample!(mpc_solver)
# We should be able to solve the MPC steps in less than mpc_solver.opts.max_δ, on average.
@test mpc_solver.stats.time <= mpc_solver.opts.iterations*mpc_solver.opts.max_δt

mpc_solver.stats.time















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
tf = 3.0  # final tim
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2
              -0.80, -0.90, # x
              -0.05, -0.30, # y
			   0.00, pi/12, # θ
			   0.60,  0.63, # v
               ]
xf = @SVector [# p1   # p2
               1.10,  0.90, # x
              -0.05, -0.05, # y
			   0.00,  0.00, # θ
			   0.60,  0.60, # v
              ]
dxf = @SVector [# p1   # p2
			   1.00,   1.00, # x
			   0.00,   0.00, # y
			   0.00,   0.00, # θ
			   0.00,   0.00, # v
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
road_length = 6.0
road_width = 0.34
ramp_length = 3.2
ramp_angle = pi/12
scenario = MergingScenario(road_length, road_width, ramp_length,
	ramp_angle, actors_radii, actors_types)

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
add_scenario_constraints(algames_conSet, scenario, px, con_inds; constraint_type=:constraint)

# Create problems
algames_ramp_merging_2_players_mpc_prob = GameProblem(model, obj, xf, tf, constraints=algames_conSet, x0=x0, N=N)

# Create solvers
algames_ramp_merging_2_players_mpc_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=0,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_ramp_merging_2_players_mpc_prob, algames_ramp_merging_2_players_mpc_opts)

# Create MPC solvers
state_noise = 5. * @SVector ([
	0.008, 0.008, 2*pi/72, 0.03, #+-50cm, +-50cm, +-25deg, +-12.5% per second
	0.008, 0.008, 2*pi/72, 0.03])
ramp_merging_2_players_mpc_opts = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=10,
	N_mpc=50,
	mpc_tf=100.0,
	min_δt=0.005,
	dxf=dxf,
	noise=state_noise)
algames_ramp_merging_2_players_mpc_solver = MPCGamesSolver(algames_solver, ramp_merging_2_players_mpc_opts)
# ilqgames_ramp_merging_2_players_mpc_solver = MPCGamesSolver(ilqgames_solver, dxf, opts_mpc)
# reset!(algames_ramp_merging_2_players_mpc_solver, reset_type=:full)
# reset!(ilqgames_ramp_merging_2_players_mpc_solver, reset_type=:full)
# solve!(algames_ramp_merging_2_players_mpc_solver; wait=false)
# resample!(algames_ramp_merging_2_players_mpc_solver)




#
#
# vis=AG.Visualizer()
# anim=AG.MeshCat.Animation()
# open(vis)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(mpc_solver, scenario;
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

















using Test

algames_prob = GameProblems.algames_ramp_merging_2_players_mpc_prob
algames_opts = GameProblems.algames_ramp_merging_2_players_mpc_opts
mpc_opts = GameProblems.ramp_merging_2_players_mpc_opts

# Create solver
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

# We solve the problem once before to avoid a substantial part of
# Julia's precompilation time for the MPC test.
solve!(algames_solver)
reset!(algames_solver, reset_type=:full)

# Create MPC solver.
mpc_solver = MPCGamesSolver(algames_solver, mpc_opts)
mpc_solver.solver.opts.log_level = Logging.Warn
mpc_solver.opts.max_δt = 0.05
reset!(mpc_solver, reset_type=:full)
mpc_solver.opts.iterations = 2
solve!(mpc_solver; wait=true)
resample!(mpc_solver)
# We should be able to solve the MPC steps in less than mpc_solver.opts.max_δ, on average.
@test mpc_solver.stats.time <= mpc_solver.opts.iterations*mpc_solver.opts.max_δt











using Test
using MeshCat

straight_prob = GameProblems.algames_straight_2_players_prob
ramp_merging_prob = GameProblems.algames_ramp_merging_2_players_prob
t_intersection_prob = GameProblems.algames_t_intersection_2_players_prob

straight_scenario = GameProblems.straight_2_players_scenario
ramp_merging_scenario = GameProblems.ramp_merging_2_players_scenario
t_intersection_scenario = GameProblems.t_intersection_2_players_scenario

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)

straight_solver = DirectGamesSolver(straight_prob, opts)
ramp_merging_solver = DirectGamesSolver(ramp_merging_prob, opts)
t_intersection_solver = DirectGamesSolver(t_intersection_prob, opts)

solve!(straight_solver)
solve!(ramp_merging_solver)
solve!(t_intersection_solver)

vis = MeshCat.Visualizer()
anim = MeshCat.Animation()
# open(vis)
# Execute this line after the MeshCat tab is open
vis, anim = animation(straight_solver, straight_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

vis, anim = animation(ramp_merging_solver, ramp_merging_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

vis, anim = animation(t_intersection_solver, t_intersection_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

@test converged(straight_solver)
@test converged(ramp_merging_solver)
@test converged(t_intersection_solver)




















using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Instantiate dynamics model
model = DoubleIntegratorGame(p=3)
n,m,pu,p = size(model)
T = Float64
px = model.px

# Discretization info
tf = 3.0  # final time
N = 21    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2   # p3
              -0.80, -1.00, -0.90, # x
              -0.05, -0.05, -0.31, # y
			  # 0.90,  0.90,  0.95, # θ
			  0.60,  0.60,  0.65, # θ
			   0.00,  0.00,  0.00, # v
               ]
xf = @SVector [# p1   # p2   # p3
               1.10,  0.70,  0.90, # x
              -0.05, -0.05, -0.05, # y
			  0.60,  0.60,  0.60, # θ
			  # 0.90,  0.90,  0.90, # θ
			   0.00,  0.00,  0.00, # v
              ]

diag_Q = [SVector{n}([1.,  0.,  0.,
					  1.,  0.,  0.,
					  1.,  0.,  0.,
					  1.,  0.,  0.]),
	      SVector{n}([0.,  1.,  0.,
		  			  0.,  1.,  0.,
					  0.,  1.,  0.,
					  0.,  1.,  0.]),
		  SVector{n}([0.,  0.,  1.,
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
inflated_actors_radii = [3.0*actor_radius for i=1:p]
actors_types = [:car for i=1:p]
road_length = 2.20
road_width = 0.42
ramp_length = 1.2
ramp_angle = pi/12
ramp_merging_3_players_penalty_scenario = MergingScenario(road_length,
	road_width, ramp_length, ramp_angle, actors_radii, actors_types)

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
con_inds = 1:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px,
	p, con_inds; constraint_type=:constraint)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, ramp_merging_3_players_penalty_scenario,
	px, con_inds; constraint_type=:constraint)

algames_ramp_merging_3_players_penalty_prob = GameProblem(model, obj, xf, tf,
	constraints=algames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
	optimality_constraint_tolerance=1e-2,
	μ_penalty=0.05,
    log_level=TO.Logging.Debug)
algames_ramp_merging_3_players_penalty_solver =
	DirectGamesSolver(
	algames_ramp_merging_3_players_penalty_prob,
	algames_opts)

# add penalty constraints
add_collision_avoidance(algames_ramp_merging_3_players_penalty_solver.penalty_constraints,
    inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)

reset!(algames_ramp_merging_3_players_penalty_solver, reset_type=:full)
@time timing_solve(algames_ramp_merging_3_players_penalty_solver)

visualize_trajectory_car(algames_ramp_merging_3_players_penalty_solver)

# using MeshCat
# vis = MeshCat.Visualizer()
# anim = MeshCat.Animation()
# open(vis)
# sleep(1.0)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(algames_ramp_merging_3_players_penalty_solver,
# 	ramp_merging_3_players_penalty_scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=false)
