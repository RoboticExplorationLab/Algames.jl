export
	MonteCarloSampler,
	MonteCarloSamplerStats,
	MonteCarloSamplerOptions,
	reset!,
	record_sample,
	num_converged,
	get_trajectory,
	get_objective,
	get_model,
	get_initial_state


@with_kw mutable struct MonteCarloSamplerStats{T,n,m,L1}
    iterations::Int = 0
    solve_time::Vector{T} = zeros(0)
    cmax::Vector{T} = zeros(0)
    H_cond::Vector{T} = zeros(0)
	iterations_total::Vector{Int} = zeros(Int,0)
	optimality_merit::Vector{T} = zeros(0)
	trajectory::Vector{Vector{KnotPoint{T,n,m,L1}}} = Vector{Vector{KnotPoint{T,n,m,L1}}}()
end

function reset!(stats::MonteCarloSamplerStats{T,n,m,L1}, S=0) where {T,n,m,L1}
    stats.iterations = 0
    stats.solve_time = zeros(S)
    stats.cmax = zeros(S)
    stats.H_cond = zeros(S)
    stats.iterations_total = zeros(Int,S)
	stats.optimality_merit = zeros(S)
	stats.trajectory = [Vector{KnotPoint{T,n,m,L1}}() for i=1:S]
    return nothing
end

@with_kw mutable struct MonteCarloSamplerOptions{n,T}
    # Options

    "Print summary at each iteration."
    verbose::Bool=false

    "Live plotting."
    live_plotting::Symbol=:off # :state, :control

    "Save samples."
    save_sampling::Bool=false

    "Number of Monte Carlo samples."
    iterations::Int = 10

	"Compute H_ condition number."
	record_condition::Bool=false

	"Amplitude of the uniformly distributed noise added to the dynamics."
    noise::SVector{n,T}=SVector{n}(zeros(n))

	log_level::Base.CoreLogging.LogLevel = TO.InnerLoop
end

struct MonteCarloSampler{T,n,m,L1}
    solver::TO.AbstractSolver{T}
    opts::MonteCarloSamplerOptions{n,T}
    stats::MonteCarloSamplerStats{T,n,m,L1}
    x0::SVector{n,T}
	xf::SVector{n,T}
    tf::T
	Q::Vector{SVector{n,T}}
	R::Vector{SVector{m,T}}
	Qf::Vector{SVector{n,T}}
    logger::TO.SolverLogger

    function MonteCarloSampler{T,n,m,L1}(
        solver::TO.AbstractSolver{T},
        opts::MonteCarloSamplerOptions{n,T},
        stats::MonteCarloSamplerStats{T,n,m,L1},
        x0::SVector{n,T},
		xf::SVector{n,T},
        tf::T,
		Q::Vector{SVector{n,T}},
		R::Vector{SVector{m,T}},
		Qf::Vector{SVector{n,T}},
        logger) where {T,n,m,L1}
        new{T,n,m,L1}(
            solver,
            opts,
            stats,
            x0,
			xf,
            tf,
			Q,
			R,
			Qf,
            logger)
    end
end

function MonteCarloSampler(solver::TO.AbstractSolver{T}, opts=MonteCarloSamplerOptions{ln,T}()) where {ln,T}
	n,m,pu,p = size(solver.model)

	# Getting the cost matrices from the solver
	# /!\ We are only getting the first and last time
	# step cost we assume that the cost are not changing with the time steps
	Q = [diag(solver.obj[i].cost[1].Q) for i=1:p]
	Ri = [zeros(m) for i=1:p]
	for i=1:p
		Ri[i][pu[i]] += diag(solver.obj[i].cost[1].R)
	end
	R = [SVector{m}(Ri[i]) for i=1:p]
	Qf = [diag(solver.obj[i].cost[end].Q) for i=1:p]


	# Init solver statistics
    stats = MonteCarloSamplerStats{T,n,m,n+m}()
    x0 = copy(solver.x0)
    xf = copy(solver.xf)
    tf = solver.tf
    logger = TO.default_logger(opts.verbose)
    sampler = MonteCarloSampler{T,n,m,n+m}(
        solver,
        opts,
        stats,
        x0,
		xf,
        tf,
		Q,
		R,
		Qf,
        logger)
    reset!(sampler)
    return sampler
end

function reset!(sampler::MonteCarloSampler{T}; reset_stats=true) where T
    if reset_stats
        reset!(sampler.stats, sampler.opts.iterations)
    end
    reset!(sampler.solver, reset_type=:full)
    return nothing
end

function record_sample(sampler::MonteCarloSampler{T},
    δt::T) where {T}
	j = sampler.solver.stats.iterations
	sampler.stats.iterations += 1
	i = sampler.stats.iterations::Int
	sampler.stats.trajectory[i] = TO.copy(sampler.solver.Z) ######################

	if string(typeof(sampler.solver).name) == "DirectGamesSolver"
		k = sampler.solver.stats.iterations_inner[j]
		cmax = sampler.solver.stats.cmax[j]
		iterations_total = sampler.solver.stats.iterations_total
		optimality_merit = sampler.solver.stats.optimality_merit[j][k]

		if sampler.opts.record_condition
			H_cond = cond(Array(sampler.solver.H_))
		else
			H_cond = NaN
		end

	    sampler.stats.solve_time[i] = δt
	    sampler.stats.cmax[i] = cmax
	    sampler.stats.H_cond[i] = H_cond
		sampler.stats.iterations_total[i] = iterations_total
		sampler.stats.optimality_merit[i] = optimality_merit
	elseif string(typeof(sampler.solver).name) == "PenaltyiLQGamesSolver"
		cmax = sampler.solver.stats.cmax[j]
		iterations_total = sampler.solver.stats.iterations
		# optimality_merit = sampler.solver.stats.optimality_merit[j][k]

	    sampler.stats.solve_time[i] = δt
	    sampler.stats.cmax[i] = cmax
		sampler.stats.iterations_total[i] = iterations_total
		# sampler.stats.optimality_merit[i] = optimality_merit
	end
    return nothing
end

function num_converged(sampler::MonteCarloSampler)
    num_converged(sampler.stats, sampler.solver)
end

function num_converged(stats::MonteCarloSamplerStats, solver::DirectGamesSolver)
    con = stats.cmax .<= solver.opts.constraint_tolerance
    opt = stats.optimality_merit .<= solver.opts.optimality_constraint_tolerance
    out = sum(con .== opt .== 1)
    return out
end

function num_converged(stats::MonteCarloSamplerStats, solver::PenaltyiLQGamesSolver)
    con = stats.cmax .<= solver.opts.constraint_tolerance
    opt = stats.optimality_merit .<= solver.opts.gradient_norm_tolerance
    out = sum(con .== opt .== 1)
    return out
end


Base.size(sampler::MonteCarloSampler{T}) where {T} = size(sampler.solver)
@inline TO.get_trajectory(sampler::MonteCarloSampler) = sampler.solver.Z
@inline TO.get_objective(sampler::MonteCarloSampler) = sampler.solver.obj
@inline TO.get_model(sampler::MonteCarloSampler) = sampler.solver.model
@inline TO.get_initial_state(sampler::MonteCarloSampler) = sampler.x0



function MonteCarloSampler(sampler_::MonteCarloSampler{T,n,m,L1},
	obj::Vector{O}, x0::SVector{n,T}, xf::SVector{n,T}) where {T,n,m,L1,O}

	if string(typeof(sampler_.solver).name) == "DirectGamesSolver"
		sampler = MonteCarloSampler{T,n,m,L1}(
			DirectGamesSolver(sampler_.solver, obj, x0, xf, sampler_.tf),
			sampler_.opts,
			sampler_.stats,
			sampler_.x0,
			sampler_.xf,
			sampler_.tf,
			sampler_.Q,
			sampler_.R,
			sampler_.Qf,
			sampler_.logger)
	elseif string(typeof(sampler_.solver).name) == "PenaltyiLQGamesSolver"
		sampler = MonteCarloSampler{T,n,m,L1}(
			PenaltyiLQGamesSolver(sampler_.solver, obj, x0, xf, sampler_.tf),
			sampler_.opts,
			sampler_.stats,
			sampler_.x0,
			sampler_.xf,
			sampler_.tf,
			sampler_.Q,
			sampler_.R,
			sampler_.Qf,
			sampler_.logger)
	end
    return sampler
end
