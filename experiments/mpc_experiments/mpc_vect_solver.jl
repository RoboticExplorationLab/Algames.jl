struct MPCVectGamesSolver11{T,S,n,m,L} <: ConstrainedSolver{T}
    # Model + Objective
	solver::Vector{S}
    opts::MPCGamesSolverOptions{n,T}
    stats::MPCGamesStats{T,n,m,L}
    x0::Vector{SVector{n,T}}
	xf::Vector{SVector{n,T}}
	dxf::SVector{n,T}
	tf::T

	Q::Vector{Vector{SVector{n,T}}}
	R::Vector{Vector{SVector{m,T}}}
	Qf::Vector{Vector{SVector{n,T}}}

    Z::Vector{KnotPoint{T,n,m,L}}
    logger::TO.SolverLogger

    function MPCVectGamesSolver11{T,S,n,m,L}(
		solver::Vector{S},
        opts::MPCGamesSolverOptions{n,T},
        stats::MPCGamesStats{T,n,m,L},
        x0::Vector{SVector{n,T}},
		xf::Vector{SVector{n,T}},
		dxf::SVector{n,T},
		tf::T,
		Q::Vector{Vector{SVector{n,T}}},
		R::Vector{Vector{SVector{m,T}}},
		Qf::Vector{Vector{SVector{n,T}}},
        Z::Vector{KnotPoint{T,n,m,L}},
        logger) where {T,S,n,m,L}
        new{T,S,n,m,L}(
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

function MPCVectGamesSolver11(solver::V, opts=MPCGamesSolverOptions{TT}()) where {V,TT}
	n,m,N = size(solver[1])
	n,m,pu,p = size(solver[1].model)

	# Getting the cost matrices from the solver
	# /!\ We are only getting the first and last time
	# step cost we assume that the cost are not changing with the time steps
	Q = [[diag(solver[i].obj[j].cost[1].Q) for j=1:p] for i=1:p]
	Ri = [[zeros(m) for j=1:p] for i=1:p]
	for i=1:p
		for j = 1:p
			Ri[i][j][pu[i]] += diag(solver[i].obj[j].cost[1].R)
		end
	end
	R = [[SVector{m}(Ri[i][j]) for j=1:p] for i=1:p]
	Qf = [[diag(solver[i].obj[j].cost[end].Q) for j=1:p] for i=1:p]

	# Init solver statistics
	S = eltype(solver)
	T = eltype(Q[1][1])
	stats = MPCGamesStats{T,n,m,n+m}()
    x0 = [copy(solver[i].x0) for i=1:p]
    xf = [copy(solver[i].xf) for i=1:p]
	dxf = copy(opts.dxf)
    tf = solver[1].tf
    L = n+m
    dt = tf/(N-1)
    Z = Traj(n,m,dt,opts.N_mpc)
    logger = TO.default_logger(opts.verbose)
	mpc_solver = MPCVectGamesSolver11{T,S,n,m,L}(
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

function AG.reset!(solver::MPCVectGamesSolver11{T}; reset_stats=true, reset_type=:nominal) where T
    if reset_stats
        reset!(solver.stats, solver.opts.iterations)
    end
	for sol in solver.solver
	   reset!(sol, reset_type=reset_type)
	end
	# reset the x0 of the solver.solver to x0 of the mpcsolver @@@@@@
    return nothing
end


function MPCVectGamesSolver11(solver_::MPCVectGamesSolver11{T,S,n,m,L},
	obj::Vector{Vector{O}},
	x0::Vector{SVector{n,T}},
	xf::Vector{SVector{n,T}},
	tf::Vector{T}=[sol.tf for sol in solver_.solver]
	) where {T,S,n,m,L,O}
    solver = MPCVectGamesSolver11{T,S,n,m,L}(
		[DirectGamesSolver(sol, obj[i], x0[i], xf[i], tf[i]) for (i,sol) in enumerate(solver_.solver)],
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

Base.size(solver::MPCVectGamesSolver11{T,S}) where {T,S} = size(solver.solver[1])
@inline TO.get_trajectory(solver::MPCVectGamesSolver11) = solver.Z
@inline TO.get_objective(solver::MPCVectGamesSolver11) = solver.solver[1].obj
@inline TO.get_model(solver::MPCVectGamesSolver11) = solver.solver[1].model
@inline TO.get_initial_state(solver::MPCVectGamesSolver11) = solver.x0
