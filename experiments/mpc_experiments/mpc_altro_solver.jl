# export
# 	MPCGamesStats,
# 	MPCGamesSolverOptions,
# 	MPCGamesSolver,
# 	reset!,
# 	get_trajectory,
# 	get_objective,
# 	get_model,
# 	get_initial_state
#
# @with_kw mutable struct MPCGamesStats{T,n,m,L}
#     iterations::Int = 0
#     solver_iterations::Vector{Int} = zeros(Int,0)
#     time::T = 0.
#     solve_time::Vector{T} = zeros(T,0)
#     cmax::Vector{T} = zeros(T,0)
# 	x0::Vector{KnotPoint{T,n,m,L}} = Vector{KnotPoint{T,n,m,L}}()
# end
#
# function reset!(stats::MPCGamesStats{T,n,m,L}, H=0) where {T,n,m,L}
#     stats.iterations = 0
#     stats.solver_iterations = zeros(Int,H)
#     stats.time = 0.
#     stats.solve_time = zeros(Int,H)
#     stats.cmax = zeros(Int,H)
# 	stats.x0 = Vector{KnotPoint{T,n,m,L}}()
# end
#
# @with_kw mutable struct MPCGamesSolverOptions{n,T} <: TO.AbstractSolverOptions{T}
#     # Options
#
#     "Print summary at each iteration."
#     verbose::Bool=false
#
#     "Live plotting."
#     live_plotting::Symbol=:off # :state, :control
#
# 	"Displacement of the goal state xf with time. xf(t) = xf(0) + t*dxf"
#     dxf::SVector{n,T}=SVector{n}(zeros(n))
#
# 	"Amplitude of the uniformly distributed noise added to the dynamics."
#     noise::SVector{n,T}=SVector{n}(zeros(n))
#
# 	"Indices of the states that are acting ignoring all others players."
# 	selfish_inds::Vector{Int}=zeros(Int,0)
#
# 	"Constant first derivative of the state ignoring all other players."
# 	selfish_dx::Vector{T}=zeros(T,0)
#
#     "Duration of simulation of the MPC."
#     mpc_tf::T=4.0
#
#     "Min duration between two MPC solves."
#     min_δt::T=0.01
#
#     "Max duration between two MPC solves, only useful for testing."
#     max_δt::T=Inf
#
# 	"Number of MPC updates to rollout the MPC."
#     N_mpc::Int=400
#
#     "Compute H_ condition number."
#     compute_cond::Bool=false
#
#     "type of game theoretic equilibrium, Nash or Stackelberg."
#     eq_type::Symbol = :nash # :nash :stackelberg
#
#     "type of information
#      pattern, Memoryless Perfect State (MPS) leads to a feedback equilibrium, Open-Loop (OL) leads to an open-loop equilibrium."
#     info_pattern::Symbol = :feedback # :feedback :open_loop
#
#     "MPCGames iterations."
#     iterations::Int = 300
#
#     "MPCGames inner iterations."
#     inner_iterations::Int = 20
#
#     log_level::Base.CoreLogging.LogLevel = TO.InnerLoop
# end

struct MPCALTROGamesSolver{T,n,m,L} <: ConstrainedSolver{T}
    # Model + Objective
    solver::Vector{TO.AbstractSolver}
    opts::MPCGamesSolverOptions{n,T}
    stats::MPCGamesStats{T,n,m,L}
    x0::SVector{n,T}
	xf::SVector{n,T}
	dxf::SVector{n,T}
	tf::T

	Q::Vector{SVector{n,T}}
	R::Vector{SVector{m,T}}
	Qf::Vector{SVector{n,T}}

    Z::Vector{KnotPoint{T,n,m,L}}
    logger::TO.SolverLogger

    function MPCALTROGamesSolver{T,n,m,L}(
        solver::Vector{TO.AbstractSolver},
        opts::MPCGamesSolverOptions{n,T},
        stats::MPCGamesStats{T,n,m,L},
        x0::SVector{n,T},
		xf::SVector{n,T},
		dxf::SVector{n,T},
		tf::T,
		Q::Vector{SVector{n,T}},
		R::Vector{SVector{m,T}},
		Qf::Vector{SVector{n,T}},
        Z::Vector{KnotPoint{T,n,m,L}},
        logger) where {T,n,m,L}
        new{T,n,m,L}(
            solver,
            opts,
            stats,
            x0,
			xf,
			dxf,
            tf,
			Q,
			R,
			Qf,
            Z,
            logger)
    end
end

function MPCALTROGamesSolver(solver::Vector{<:TO.AbstractSolver{T}}, opts=MPCGamesSolverOptions{T}()) where {T}
	n,m,N = size(solver[1])
	model = get_model(solver[1])
	p = model.p
	pu = model.pu

	# Getting the cost matrices from the solver
	# /!\ We are only getting the first and last time
	# step cost we assume that the cost are not changing with the time steps
	Q = [diag(TO.get_objective(solver[i]).obj.cost[1].Q) for i=1:p]
	Ri = [zeros(m) for i=1:p]
	for i=1:p
		Ri[i][pu[i]] += diag(TO.get_objective(solver[i]).obj.cost[1].R[pu[i],pu[i]])
	end
	R = [SVector{m}(Ri[i]) for i=1:p]
	Qf = [diag(TO.get_objective(solver[i]).obj.cost[end].Q) for i=1:p]


	# Init solver statistics
    stats = MPCGamesStats{T,n,m,n+m}()
    x0 = copy(TO.get_initial_state(solver[1].solver_al))
    xf = copy(solver[1].solver_al.solver_uncon.xf)
	dxf = copy(opts.dxf)
    tf = solver[1].solver_al.solver_uncon.tf
    L = n+m
    dt = tf/(N-1)
    Z = Traj(n,m,dt,opts.N_mpc)
    logger = TO.default_logger(opts.verbose)
    mpc_solver = MPCALTROGamesSolver{T,n,m,L}(
        solver,
        opts,
        stats,
        x0,
		xf,
		dxf,
        tf,
		Q,
		R,
		Qf,
        Z,
        logger)
    reset!(mpc_solver)
    return mpc_solver
end

function AG.reset!(solver::MPCALTROGamesSolver{T}; reset_stats=true, reset_type=:nominal) where T
    if reset_stats
        reset!(solver.stats, solver.opts.iterations)
    end
	for sol in solver.solver
		reset!(sol, reset_type=reset_type)
	end
	# reset the x0 of the solver.solver to x0 of the mpcsolver @@@@@@
    return nothing
end

function MPCALTROGamesSolver(solver_::MPCALTROGamesSolver{T,n,m,L},
	obj::Vector{O},
	x0::SVector{n,T},
	xf::SVector{n,T},
	tf::T=solver_.solver.tf
	) where {T,n,m,L,O}
    solver = MPCALTROGamesSolver{T,n,m,L}(
		DirectGamesSolver(solver_.solver, obj, x0, xf, tf),
		solver_.opts,
		solver_.stats,
		solver_.x0,
		solver_.xf,
		solver_.dxf,
		solver_.tf,
		solver_.Q,
		solver_.R,
		solver_.Qf,
		solver_.Z,
		solver_.logger)
    return solver
end

Base.size(solver::MPCALTROGamesSolver{T}) where {T} = size(solver.solver[1])
@inline TO.get_trajectory(solver::MPCALTROGamesSolver) = solver.Z
@inline TO.get_objective(solver::MPCALTROGamesSolver) = [get_objective(sol.obj) for sol in solver.solver]
@inline TO.get_model(solver::MPCALTROGamesSolver) = solver.solver[1].model
@inline TO.get_initial_state(solver::MPCALTROGamesSolver) = solver.x0
