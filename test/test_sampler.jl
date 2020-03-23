using Test

algames_solver = GameProblems.algames_ramp_merging_3_players_solver
algames_solver.opts.min_steps_per_iteration = 1
ilqgames_solver = GameProblems.ilqgames_ramp_merging_3_players_solver

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
# algames_sampler.stats.cmax
# ilqgames_sampler.stats
# mean(algames_sampler.stats.solve_time)
# mean(ilqgames_sampler.stats.solve_time)
#
# algames_sampler.stats.cmax
# algames_sampler.stats.optimality_merit
#
# all(algames_sampler.stats.cmax .<= algames_solver.opts.constraint_tolerance)
# all(algames_sampler.stats.optimality_merit .<= algames_solver.opts.optimality_constraint_tolerance)
#
# all(ilqgames_sampler.stats.cmax .<= ilqgames_solver.opts.constraint_tolerance)
# all(ilqgames_sampler.stats.optimality_merit .<= ilqgames_solver.opts.gradient_norm_tolerance)
