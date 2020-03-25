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
mpc_solver.opts.iterations = 10
solve!(mpc_solver; wait=true)
resample!(mpc_solver)
# We should be able to solve the MPC steps in less than mpc_solver.opts.max_δ, on average.
@test mpc_solver.stats.time <= mpc_solver.opts.iterations*mpc_solver.opts.max_δt
