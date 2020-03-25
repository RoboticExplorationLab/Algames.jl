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
