using Test

ramp_merging_2_players_prob = GameProblems.algames_ramp_merging_2_players_prob
# ramp_merging_3_players_prob = GameProblems.algames_ramp_merging_3_players_prob
# ramp_merging_4_players_prob = GameProblems.algames_ramp_merging_4_players_prob
straight_2_players_prob = GameProblems.algames_straight_2_players_prob
t_intersection_2_players_prob = GameProblems.algames_t_intersection_2_players_prob
# t_intersection_3_players_prob = GameProblems.algames_t_intersection_3_players_prob
# t_intersection_4_players_prob = GameProblems.algames_t_intersection_4_players_prob

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)

ramp_merging_2_players_solver = DirectGamesSolver(ramp_merging_2_players_prob, opts)
# ramp_merging_3_players_solver = DirectGamesSolver(ramp_merging_3_players_prob, opts)
# ramp_merging_4_players_solver = DirectGamesSolver(ramp_merging_4_players_prob, opts)
straight_2_players_solver = DirectGamesSolver(straight_2_players_prob, opts)
t_intersection_2_players_solver = DirectGamesSolver(t_intersection_2_players_prob, opts)
# t_intersection_3_players_solver = DirectGamesSolver(t_intersection_3_players_prob, opts)
# t_intersection_4_players_solver = DirectGamesSolver(t_intersection_4_players_prob, opts)

# Solve with ALGAMES
solve!(ramp_merging_2_players_solver)
# solve!(ramp_merging_3_players_solver)
# solve!(ramp_merging_4_players_solver)
solve!(straight_2_players_solver)
solve!(t_intersection_2_players_solver)
# solve!(t_intersection_3_players_solver)
# solve!(t_intersection_4_players_solver)

@test converged(ramp_merging_2_players_solver)
# @test converged(ramp_merging_3_players_solver)
# @test converged(ramp_merging_4_players_solver)
@test converged(straight_2_players_solver)
@test converged(t_intersection_2_players_solver)
# @test converged(t_intersection_3_players_solver)
# @test converged(t_intersection_4_players_solver)

# # Check the gradient and Hessian
# solver_2_players.stats.α
# solver_3_players.stats.α
#
# visualize_trajectory_car(solver_2_players)
