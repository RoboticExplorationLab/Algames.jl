export
	solve!,
	update_traj!,
	step!,
	need_resolve,
	generate_new_x0,
	mpc_propagate_dynamics,
	record_iteration!,
	evaluate_convergence,
	resample!


# Generic solve methods
"MPCGames solve method (non-allocating)"
function TO.solve!(solver::MPCVectGamesSolver11{T}; wait::Bool=false, fix_seed::Bool=true) where {T<:AbstractFloat}
	fix_seed ? Random.seed!(100) : nothing
    TO.set_verbosity!(solver.opts)
    TO.clear_cache!(solver.opts)
    solver.stats.iterations = 0
	n,m,N = size(solver)
	n,m,pu,p = size(get_model(solver))

    for i = 1:solver.opts.iterations
		# TO.Logging.@info "MPC Solver iteration = ", i
        dt_step = @elapsed δt, new_x0, new_xf = AG.step!(solver)
		if wait
			sleep(0.2)
		end
		# Update the solver.
		dt_new_obj = @elapsed new_obj = [[TO.LQRObjective(
								Diagonal(solver.Q[i][j]),
								Diagonal(SVector{length(pu[j])}(solver.R[i][j][pu[j]])),
								Diagonal(solver.Qf[i][j]),
								new_xf[i],
								N,
								checks=false) for j=1:p] for i=1:p]
		dt_new_solver = @elapsed solver = MPCVectGamesSolver11(solver, new_obj, new_x0, new_xf) ###############################
		# for i = 1:p
		# 	dt_update_traj = @elapsed AG.update_traj!(solver, i, δt, new_x0) ###############################################################
		# end
        dt_eval_cv = @elapsed AG.evaluate_convergence(solver) ? break : nothing
		# @show dt_step dt_new_obj dt_new_solver dt_update_traj dt_eval_cv
    end
    return solver
end

function AG.update_traj!(solver::MPCVectGamesSolver11{T}, player_id::Int,  δt::T, x0::SVector{ln,T}) where {ln,T}
	#need to set the correct x0
	# need to update the controls in Z
	# don't need to update Zbar
	# Set the initial state and controls
	n,m,N = size(solver)
	Z = solver.solver[player_id].Z
	t = δt
	dt = 0.
	k = 0
	while t > 0 && k < N-1
		k +=1
		dt = min(Z[k].dt, t)
		t -= dt
	end
	ulast = TO.control(Z[k])
	unext = TO.control(Z[k+1])
	# Linear interpolation between last applied control and the next one
	u0 = ulast + (unext - ulast)*dt/Z[k].dt
	# Update first knotpoint
	Z[1] = TO.KnotPoint(x0, u0, Z[k].dt-dt, 0.)
	# Update the N-k-1 middle knotpoints
	for (j,l) in enumerate(2:N-k)
		x = TO.state(Z[k+j])
		u = TO.control(Z[k+j])
		dt = Z[k+j].dt
		Z[l] = TO.KnotPoint(x,u,dt)
	end
	# Update the k-1 final knotpoints
	for l = N-k+1:N-1
		# Does not really matter since the solver rolls the traj out.
		x = TO.state(Z[N-1])
		# We pad the last part of the control trajectory
		# with the last control of the previous trajectory.
		u = TO.control(Z[N-1])
		dt = Z[N-1].dt
		Z[l] = TO.KnotPoint(x,u,dt)
	end
	# Update the final knotpoint
	Z[N] = TO.KnotPoint(TO.state(Z[N]),m)
	return nothing
end


"""
Take one step of MPCGames algorithm (non-allocating)
"""
function AG.step!(solver::MPCVectGamesSolver11{T}) where {T}
	# δt = @elapsed resolve = need_resolve(solver.solver)
	δt = @elapsed resolve = true
	if resolve
	    δt += @elapsed solve!.(solver.solver) ####################
		TO.Logging.@info δt
		δt = clamp(δt, solver.opts.min_δt, solver.opts.max_δt)
		record_iteration!(solver, δt; resolve=true)
	# else
	# 	TO.Logging.@info δt
	# 	δt = clamp(δt, solver.opts.min_δt, solver.opts.max_δt)
	# 	record_iteration!(solver, δt; resolve=false)
	end

	# Add noise and propagate dynamics forward
	new_x0 = generate_new_x0(solver, δt)
	new_x0 = [new_x0 for i=1:p]

	time_remaining = solver.opts.mpc_tf - solver.stats.time - δt
	new_xf = [sol.xf + solver.dxf*δt for sol in solver.solver]

	TO.Logging.@info time_remaining
	if solver.opts.live_plotting == :on
		for sol in solver.solver
			visualize_trajectory_car(sol;save_figure=false)
		end
	end
	for sol in solver.solver
		AG.reset!(sol, reset_type=:mpc)
	end
    return δt, new_x0, new_xf
end

function AG.generate_new_x0(solver::MPCVectGamesSolver11, δt::T) where T
	# Add noise and propagate dynamics forward

	n,m,N = size(solver)
	n,m,pu,p = size(get_model(solver))
	x0s = [mpc_propagate_dynamics(sol.model, sol.Z, δt) for sol in solver.solver]
	new_x0 = zeros(n)
	pxs = [[(j-1)*p+i for j=1:Int(n/p)] for i=1:p]
	for i = 1:p
		new_x0[pxs[i]] = x0s[i][pxs[i]]
	end

	selfish_x0 = Array(copy(new_x0))
	for (i, ind) in enumerate(solver.opts.selfish_inds)
		selfish_x0[ind] += δt*solver.opts.selfish_dx[i]
	end
	rnd = SVector{n}(rand(n) .- 0.5)
	new_x0 += solver.opts.noise .* rnd .* δt
	for i in setdiff(1:n, solver.opts.selfish_inds)
		selfish_x0[i] = new_x0[i]
	end
	new_x0 = SVector{n}(selfish_x0)
	return new_x0
end

"""
Stash iteration statistics
"""
function AG.record_iteration!(solver::MPCVectGamesSolver11, δt; resolve::Bool=true)
    solver.stats.iterations += 1
    i = solver.stats.iterations::Int
	solver.stats.time += δt
	solver.stats.solve_time[i] = δt

	n,m,N = size(solver)
	n,m,pu,p = size(get_model(solver))

	z0 = merge_z0(solver.solver, [sol.Z for sol in solver.solver], δt)
	push!(solver.stats.x0, z0)

	if resolve
		solver.stats.solver_iterations[i] = solver.solver[1].stats.iterations_total
		j = solver.solver[1].stats.iterations
		solver.stats.cmax[i] = maximum([solver.solver[i].stats.cmax[j] for i=1:p])
	# else
	# 	solver.stats.solver_iterations[i] = 0
	# 	TO.evaluate!(solver.solver.constraints, solver.solver.Z)
	# 	TO.evaluate!(solver.solver.dyn_constraints, solver.solver.Z)
	# 	TO.max_violation!(solver.solver.constraints)
	# 	TO.max_violation!(solver.solver.dyn_constraints)
	# 	cmax = max(maximum(solver.solver.constraints.c_max),
	# 		maximum(solver.solver.dyn_constraints.c_max))
	# 	solver.stats.cmax[i] = cmax
	end

	cmax_collision, cmax_boundary, cmax_bound =	AG.evaluate_constraint_violation(
		solver.solver, [sol.Z for sol in solver.solver])
	solver.stats.cmax_collision[i] = cmax_collision
	solver.stats.cmax_boundary[i] = cmax_boundary
	solver.stats.cmax_bound[i] = cmax_bound

	@logmsg TO.InnerLoop :iter value=i
	@logmsg TO.InnerLoop :inner_iter value=solver.stats.solver_iterations[i]
	@logmsg TO.InnerLoop :time value=solver.stats.time
	@logmsg TO.InnerLoop :δt value=solver.stats.time[i]
    @logmsg TO.InnerLoop :cmax value=solver.stats.cmax[i]
    if solver.opts.verbose
        print_level(InnerLoop)
    end
    return nothing
end


# """
# $(SIGNATURES)
# Check convergence conditions for MPCR
# """
function AG.evaluate_convergence(solver::MPCVectGamesSolver11)
    # Get current iterations
	i = solver.stats.iterations
	t = solver.stats.time
	# Check constraint violation
	if solver.stats.cmax_collision[i] >= 2e-2
		solver.stats.failure_status = true
		TO.Logging.@info "Collision constraint is violated cmax_collision >= 2e-3."
		return true
	end
	if solver.stats.cmax_boundary[i] >= 2e-2
		solver.stats.failure_status = true
		TO.Logging.@info "Boundary constraint is violated cmax_boundary >= 2e-3."
		return true
	end
	if solver.stats.cmax_bound[i] >= 1e-1
		solver.stats.failure_status = true
		TO.Logging.@info "Bound constraint is violated cmax_bound >= 1e-2."
		return true
	end
	# Check MPC time
    if t >= solver.opts.mpc_tf
		TO.Logging.@info "MPC solver reached final time."
        return true
    end
    # Check total iterations
    if i >= solver.opts.iterations
		TO.Logging.@info "MPC solver max # iterations."
        return true
    end
    return false
end

function AG.resample!(solver::MPCVectGamesSolver11{T}) where T
	N_mpc = solver.opts.N_mpc
	dt_mpc = solver.opts.mpc_tf /(N_mpc-1)
	M = solver.stats.iterations
	ct = 1
	x0 = solver.stats.x0[1]
	solver.Z[1] = TO.KnotPoint(TO.state(x0),TO.control(x0),dt_mpc,0.)
	rest = 0.
	for j = 2:M-1
		xprev = TO.state(solver.stats.x0[j-1])
		uprev = TO.control(solver.stats.x0[j-1])
		xnext = TO.state(solver.stats.x0[j])
		unext = TO.control(solver.stats.x0[j])
		Δt = copy(solver.stats.x0[j-1].dt) + rest
		α = 0.
		while Δt > dt_mpc
			Δt -= dt_mpc
			α += (dt_mpc-rest)/copy(solver.stats.x0[j-1].dt)
			x = xprev + (xnext - xprev)*α
			u = uprev + (unext - uprev)*α
			rest = 0.
			ct += 1
			if ct > N_mpc
				break
			end
			solver.Z[ct] = TO.KnotPoint(x,u,dt_mpc,0.)
		end
		rest = Δt
	end
	return nothing
end



function AG.evaluate_constraint_violation(solvers::Vector{S}, Zs::Vector{Tr}) where {S,Tr}
	# Collision, Boundary, Bound
	cmax_collision = 0.
	cmax_boundary = 0.
	cmax_bound = 0.

	n,m,N = size(solvers[1])
	n,m,pu,p = size(get_model(solvers[1]))
	pxs = [[(j-1)*p+i for j=1:Int(n/p)] for i=1:p]


	z0 = merge_z0(solvers, Zs)

	for i = 1:p
		for con in solvers[i].constraints
			evaluate!(con.vals, con.con, [z0], [1])
			if typeof(con.con) <: CollisionConstraint
				cmax_collision = max(cmax_collision, maximum(con.vals[1]))
			elseif typeof(con.con) <: BoundaryConstraint
				cmax_boundary = max(cmax_boundary, maximum(con.vals[1]))
			elseif typeof(con.con) <: BoundConstraint
				inds = [pu[i]; m .+ pu[i]]
				cmax_bound = max(cmax_bound, maximum(con.vals[1][inds]))
			end
		end
	end
	return cmax_collision, cmax_boundary, cmax_bound
end

function merge_z0(solvers::Vector{S},
	Zs::Vector{Tr}=[sol.Z for sol in solvers], δt::T=0.0) where {S,Tr,T}

	n,m,N = size(solvers[1])
	n,m,pu,p = size(get_model(solvers[1]))
	x0s = [TO.state(solvers[i].Z[1]) for i=1:p]
	u0s = [TO.control(solvers[i].Z[1]) for i=1:p]

	x0 = zeros(n)
	u0 = zeros(m)
	pxs = [[(j-1)*p+i for j=1:Int(n/p)] for i=1:p]
	for i = 1:p
		x0[pxs[i]] = x0s[i][pxs[i]]
		u0[pu[i]] = u0s[i][pu[i]]
	end

	z0 = TO.KnotPoint(x0, u0, δt, 0.)
	return z0
end

function initial_feasibility_check(solver::DirectGamesSolver; display::Bool=false)
	# Check that the initial state is feasible.
	reset!(solver, reset_type=:full)
	TO.rollout!(solver)
	display ? visualize_trajectory_car(solver) : nothing
	TO.evaluate!(solver.constraints, solver.Z)
	cmax = 0.0
	for con in solver.constraints
		cmax = max(cmax, maximum(con.vals[1]))
	end
	feasible = cmax <= 0.0
	return feasible
end

function monte_carlo_analysis(solvers::Vector{S}, opts_mpc::MPCGamesSolverOptions{n,T},
	state_noise::SVector{n,T}, iterations::Int=2; fix_seed::Bool=true) where {S,n,T}
	fix_seed ? Random.seed!(100) : nothing

	p = get_model(solvers[1]).p
	# Copy original x0
	nominal_x0s = [deepcopy(sol.x0) for sol in solvers]
	x_pos_rank = zeros(Int, iterations)
	failure = falses(iterations)

	p = length(solvers)
	s = iterations
	while s >= 1
        # Set x0 and xf with noise
		δx0 = (SVector{n}(rand(n)) .- 0.5) .* state_noise
		noisy_x0s = [nominal_x0s[i] + δx0 for i=1:p]
		solvers = [DirectGamesSolver(sol, sol.obj, noisy_x0s[i], sol.xf) for (i,sol) in enumerate(solvers)]

		# Check that the initial state is feasible.
		feasible = all(initial_feasibility_check.(solvers, display=false))
		if !feasible
			@show "infeasible"
			continue
		end
		reset!.(solvers, reset_type=:full)

		mpc_solver = MPCVectGamesSolver11(solvers, opts_mpc)
		reset!(mpc_solver, reset_type=:full)
		solve!(mpc_solver; wait=false, fix_seed=false)
		resample!(mpc_solver)

		visualize_state(mpc_solver)
		failure[s] = mpc_solver.stats.failure_status
		if !failure[s]
			initial_state = TO.state(mpc_solver.Z[1])
			final_state = TO.state(mpc_solver.Z[end-10])
			y_pos_initial = initial_state[Int(n/p)+1:2*Int(n/p)]
			x_pos_final = final_state[1:Int(n/p)]
			id_merging = findmin(y_pos_initial)[2]
			@show id_merging
			@show y_pos_initial
			@show x_pos_final
			x_pos_rank[s] = findfirst(x->x==x_pos_final[id_merging], sort(x_pos_final, rev=true))
		end
		s -= 1
    end
	return failure, x_pos_rank
end

function monte_carlo_analysis(solver::S, opts_mpc::MPCGamesSolverOptions{n,T},
	state_noise::SVector{n,T}, iterations::Int=2; fix_seed::Bool=true) where {S,n,T}
	fix_seed ? Random.seed!(100) : nothing

	p = get_model(solver).p
	# Copy original x0
	nominal_x0 = deepcopy(solver.x0)
	x_pos_rank = zeros(Int, iterations)
	failure = falses(iterations)

	s = iterations
	while s >= 1
        # Set x0 and xf with noise
		δx0 = (SVector{n}(rand(n)) .- 0.5) .* state_noise
		noisy_x0 = nominal_x0 + δx0
		solver = DirectGamesSolver(solver, solver.obj, noisy_x0, solver.xf)

		# Check that the initial state is feasible.
		feasible = initial_feasibility_check(solver, display=false)
		if !feasible
			@show "infeasible"
			continue
		end
		reset!(solver, reset_type=:full)
		mpc_solver = MPCGamesSolver(solver, opts_mpc)
		reset!(mpc_solver, reset_type=:full)
		solve!(mpc_solver; wait=false, fix_seed=false)
		resample!(mpc_solver)

		visualize_state(mpc_solver)
		failure[s] = mpc_solver.stats.failure_status

		if !failure[s]
			initial_state = TO.state(mpc_solver.Z[1])
			final_state = TO.state(mpc_solver.Z[end-10])
			y_pos_initial = initial_state[Int(n/p)+1:2*Int(n/p)]
			x_pos_final = final_state[1:Int(n/p)]
			id_merging = findmin(y_pos_initial)[2]
			@show id_merging
			@show y_pos_initial
			@show x_pos_final
			x_pos_rank[s] = findfirst(x->x==x_pos_final[id_merging], sort(x_pos_final, rev=true))
		end
		s -= 1
    end
	return failure, x_pos_rank
end
