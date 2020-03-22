using Test

prob = copy(GameProblems.ilqgames_ramp_merging)

# Solve with iLQGames
T = Float64
ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05,
    log_level=TO.Logging.Warn,
    )
ilqgames = PenaltyiLQGamesSolver(prob, ilqgames_opts)
pen = ones(length(ilqgames.constraints))*100.0
set_penalty!(ilqgames, pen);
solve!(ilqgames)
@test converged(ilqgames)
