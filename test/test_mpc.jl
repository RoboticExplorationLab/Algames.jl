using Test


@test 1 == 1
#
# using Logging
# algames_mpc_solver = GameProblems.algames_ramp_merging_3_players_mpc_solver
# reset!(algames_mpc_solver, reset_type=:full)
# algames_mpc_solver.solver.opts.constraint_tolerance = 1e-3
# algames_mpc_solver.solver.opts.optimality_constraint_tolerance = 1e-2
# algames_mpc_solver.solver.opts.log_level = Logging.Debug
# solve!(algames_mpc_solver; wait=false)
# resample!(algames_mpc_solver)
#
#
# algames_mpc_solver.stats
# @test
#
# visualize_trajectory_car(algames_mpc_solver.solver)
# a =10
# a =10
# a =10
# a =10
# a =10



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
# mpc_solver.opts.log_level = AG.Logging.Debug
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
