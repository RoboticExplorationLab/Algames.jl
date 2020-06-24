using Statistics
using Random
using Plots
using Logging

include("algames_solver_2_players.jl")
include("approx_solver_2_players.jl")
include("const_vel_solver_2_players.jl")
include("mpc_vect_solver.jl")
include("mpc_vect_methods.jl")

# Algames Solver
ramp_merging_2_players_unicycle_penalty_scenario
algames_solver
approx_solvers
const_vel_solvers

# Define initial and final states (be sure to use Static Vectors!)
dxf = @SVector [
               0.60, 0.70,
			   0.00, 0.00,
			   0.00, 0.00,
			   0.00, 0.00,
			   ]
state_noise = 2. * SVector{n}([
	0.008,   0.008,
	0.008,   0.008,
	2*pi/72, 2*pi/72,
	0.03,    0.03]) #+-50cm, +-50cm, +-25deg, +-12.5% per second
opts_mpc = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=2000,
	N_mpc=50,
	mpc_tf=8.0,
	min_δt=0.005,
	max_δt=0.40,
    dxf=dxf,
	noise=state_noise)

algames_mpc_solver = MPCGamesSolver(algames_solver, opts_mpc)
reset!(algames_mpc_solver, reset_type=:full)
solve!(algames_mpc_solver; wait=false)
resample!(algames_mpc_solver)
algames_mpc_solver.stats.failure_status

approx_mpc_solver = MPCVectGamesSolver11(approx_solvers, opts_mpc)
reset!(approx_mpc_solver, reset_type=:full)
solve!(approx_mpc_solver; wait=false)
resample!(approx_mpc_solver)
approx_mpc_solver.stats.failure_status

const_vel_mpc_solver = MPCVectGamesSolver11(const_vel_solvers, opts_mpc)
reset!(const_vel_mpc_solver, reset_type=:full)
solve!(const_vel_mpc_solver; wait=false)
resample!(const_vel_mpc_solver)
const_vel_mpc_solver.stats.failure_status

maximum(algames_mpc_solver.stats.cmax_collision)
maximum(algames_mpc_solver.stats.cmax_boundary)
maximum(algames_mpc_solver.stats.cmax_bound)

maximum(approx_mpc_solver.stats.cmax_collision)
maximum(approx_mpc_solver.stats.cmax_boundary)
maximum(approx_mpc_solver.stats.cmax_bound)

maximum(const_vel_mpc_solver.stats.cmax_collision)
maximum(const_vel_mpc_solver.stats.cmax_boundary)
maximum(const_vel_mpc_solver.stats.cmax_bound)



mpc_state_noise = SVector{n}([
	0.20,    0.20,  # 20m
	0.02,    0.02,  # 40cm
	pi/30,   pi/30, # 6 deg
	0.05,    0.05]) # 1m/s
iter = 100
algames_failure = monte_carlo_analysis(
	algames_solver,
	opts_mpc,
	mpc_state_noise,
	iter)

approx_failure = monte_carlo_analysis(
	approx_solvers,
	opts_mpc,
	mpc_state_noise,
	iter)

const_vel_failure = monte_carlo_analysis(
	const_vel_solvers,
	opts_mpc,
	mpc_state_noise,
	iter)


a = 10
a = 10
a = 10
a = 10


sum(algames_failure)
sum(approx_failure)
sum(const_vel_failure)

5
6
8


vis=AG.Visualizer()
anim=AG.MeshCat.Animation()
open(vis)
# Execute this line after the MeshCat tab is open
vis, anim = animation(approx_mpc_solver, ramp_merging_2_players_unicycle_penalty_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)


vis, anim = animation(const_vel_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

#
# #
# # vis=AG.Visualizer()
# # anim=AG.MeshCat.Animation()
# # open(vis)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(const_vel_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
# vis=vis, anim=anim,
# open_vis=false,
# display_actors=true,
# display_trajectory=true)


# iter = mpc_solver.stats.iterations
# mean_solve_time = mean(mpc_solver.stats.solve_time[1:iter])
# update_freq = mpc_solver.stats.iterations/mpc_solver.stats.time
# std_solve_time = std(mpc_solver.stats.solve_time)
# largest_solve_time = maximum(mpc_solver.stats.solve_time)
#
# iter = const_vel_mpc_solver.stats.iterations
# mean_solve_time = mean(const_vel_mpc_solver.stats.solve_time[1:iter])
# update_freq = const_vel_mpc_solver.stats.iterations/mpc_solver.stats.time
# std_solve_time = std(const_vel_mpc_solver.stats.solve_time)
# largest_solve_time = maximum(const_vel_mpc_solver.stats.solve_time)
