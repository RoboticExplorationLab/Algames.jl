using Test

ramp_merging_2_players_prob = GameProblems.ilqgames_ramp_merging_2_players_prob
# ramp_merging_3_players_prob = GameProblems.ilqgames_ramp_merging_3_players_prob
# ramp_merging_4_players_prob = GameProblems.ilqgames_ramp_merging_4_players_prob
straight_2_players_prob = GameProblems.ilqgames_straight_2_players_prob
t_intersection_2_players_prob = GameProblems.ilqgames_t_intersection_2_players_prob
# t_intersection_3_players_prob = GameProblems.ilqgames_t_intersection_3_players_prob
# t_intersection_4_players_prob = GameProblems.ilqgames_t_intersection_4_players_prob

opts = PenaltyiLQGamesSolverOptions(
    iterations=600,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.02,
    log_level=TO.Logging.Warn,
    )

ramp_merging_pen_2 = ones(length(ramp_merging_2_players_prob.constraints))*10000.0
# ramp_merging_pen_3 = ones(length(ramp_merging_3_players_prob.constraints))*1000.0
# ramp_merging_pen_4 = ones(length(ramp_merging_4_players_prob.constraints))*1000.0
straight_pen_2 = ones(length(straight_2_players_prob.constraints))*10000.0
t_intersecion_pen_2 = ones(length(t_intersection_2_players_prob.constraints))*10000.0
# t_intersecion_pen_3 = ones(length(t_intersection_3_players_prob.constraints))*1000.0
# t_intersecion_pen_4 = ones(length(t_intersection_4_players_prob.constraints))*1000.0

ramp_merging_2_players_solver = PenaltyiLQGamesSolver(ramp_merging_2_players_prob, opts)
# ramp_merging_3_players_solver = PenaltyiLQGamesSolver(ramp_merging_3_players_prob, opts)
# ramp_merging_4_players_solver = PenaltyiLQGamesSolver(ramp_merging_4_players_prob, opts)
straight_2_players_solver = PenaltyiLQGamesSolver(straight_2_players_prob, opts)
t_intersection_2_players_solver = PenaltyiLQGamesSolver(t_intersection_2_players_prob, opts)
# t_intersection_3_players_solver = PenaltyiLQGamesSolver(t_intersection_3_players_prob, opts)
# t_intersection_4_players_solver = PenaltyiLQGamesSolver(t_intersection_4_players_prob, opts)

set_penalty!(ramp_merging_2_players_solver, ramp_merging_pen_2)
# set_penalty!(ramp_merging_3_players_solver, ramp_merging_pen_3)
# set_penalty!(ramp_merging_4_players_solver, ramp_merging_pen_4)
set_penalty!(straight_2_players_solver, straight_pen_2)
set_penalty!(t_intersection_2_players_solver, t_intersecion_pen_2)
# set_penalty!(t_intersection_3_players_solver, t_intersecion_pen_3)
# set_penalty!(t_intersection_4_players_solver, t_intersecion_pen_4)


# Solve with iLQGames
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
