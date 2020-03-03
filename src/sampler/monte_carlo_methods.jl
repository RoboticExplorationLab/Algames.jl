export
	monte_carlo_sampling


function monte_carlo_sampling(sampler::MonteCarloSampler{T}; display_threshold::Int=50) where {n1,T}
	Random.seed!(100)
	n,m,N = size(sampler)
	n,m,pu,p = size(sampler.solver.model)
	s = sampler.opts.iterations
    while s >= 1
        # Set x0 and xf with noise
		δx0 = (SVector{n}(rand(n)) .- 0.5) .* sampler.opts.noise
        noisy_x0 = sampler.x0 + δx0
		noisy_xf = sampler.xf
		noisy_obj = [TO.LQRObjective(
							Diagonal(sampler.Q[i]),
							Diagonal(SVector{length(pu[i])}(sampler.R[i][pu[i]])),
							Diagonal(sampler.Qf[i]),
							noisy_xf,
							N) for i=1:p]

        sampler = MonteCarloSampler(sampler, noisy_obj, noisy_x0, noisy_xf)

		# Check that the initial state is feasible.
		reset!(sampler.solver, reset_type=:full)
		rollout!(sampler.solver)
		TO.evaluate!(sampler.solver.constraints, sampler.solver.Z)
		cmax = 0.0
		for con in sampler.solver.constraints
			cmax = max(cmax, maximum(con.vals[1]))
		end
		if cmax > 0.0
			continue
		end

		reset!(sampler.solver, reset_type=:full)
		δt = @elapsed solve!(sampler.solver)
		# if display || sampler.solver.stats.iterations_total > display_threshold
		# if display || sampler.solver.stats.iterations > display_threshold
		if sampler.opts.live_plotting == :on
			visualize_trajectory_car(sampler.solver)
			# animation(sampler.solver, scenario)
		end
        record_sample(sampler, δt)
		s -= 1
    end
	if sampler.opts.save_sampling
		save("files/sampling/sampler_stats_" *
	    	Dates.format(now(), "HH:MM:SS:sss") *
	    	".jld2",
	    	"sampler.stats",
			sampler.stats,
			)
	end
    return nothing
end
