using Test

algames_prob = GameProblems.algames_ramp_merging_2_players_prob
opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_prob, opts)

ilqgames_prob = GameProblems.ilqgames_ramp_merging_2_players_prob
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

num_samples = 5
state_noise = @SVector [ # Uniform noise around x0
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
