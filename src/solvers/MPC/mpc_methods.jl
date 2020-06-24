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
function TO.solve!(solver::MPCGamesSolver{T}; wait::Bool=false, fix_seed::Bool=true) where {T<:AbstractFloat}
	fix_seed ? Random.seed!(100) : nothing
    TO.set_verbosity!(solver.opts)
    TO.clear_cache!(solver.opts)
    solver.stats.iterations = 0
	n,m,N = size(solver.solver)
	n,m,pu,p = size(solver.solver.model)

    for i = 1:solver.opts.iterations
		# TO.Logging.@info "MPC Solver iteration = ", i
		# println("MPC Solver iteration = ", i)
        dt_step = @elapsed δt, new_x0, new_xf = step!(solver)
		if wait
			sleep(0.2)
		end
		# Update the solver.
		dt_new_obj = @elapsed new_obj = [TO.LQRObjective(
								Diagonal(solver.Q[i]),
								Diagonal(SVector{length(pu[i])}(solver.R[i][pu[i]])),
								Diagonal(solver.Qf[i]),
								new_xf,
								N,
								checks=false) for i=1:p]

		# dt_new_solver = @elapsed solver = MPCGamesSolver(solver, solver.solver.obj, new_x0, new_xf) #######################
		dt_new_solver = @elapsed solver = MPCGamesSolver(solver, new_obj, new_x0, new_xf) ###############################
		dt_update_traj = @elapsed update_traj!(solver, δt, new_x0)
        dt_eval_cv = @elapsed evaluate_convergence(solver) ? break : nothing
		# @show dt_step dt_new_obj dt_new_solver dt_update_traj dt_eval_cv
    end
    return solver
end

function update_traj!(solver::MPCGamesSolver{T}, δt::T, x0::SVector{ln,T}) where {ln,T}
	#need to set the correct x0
	# need to update the controls in Z
	# don't need to update Zbar
	# Set the initial state and controls
	n,m,N = size(solver)
	Z = solver.solver.Z
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
function step!(solver::MPCGamesSolver{T}) where {T}
	δt = @elapsed resolve = need_resolve(solver.solver)
	if resolve
	    δt += @elapsed solve!(solver.solver) ####################
		TO.Logging.@info δt
		δt = clamp(δt, solver.opts.min_δt, solver.opts.max_δt)
		record_iteration!(solver, δt; resolve=true)
	else
		TO.Logging.@info δt
		δt = clamp(δt, solver.opts.min_δt, solver.opts.max_δt)
		record_iteration!(solver, δt; resolve=false)
	end

	# Add noise and propagate dynamics forward
	new_x0 = generate_new_x0(solver, δt)

	time_remaining = solver.opts.mpc_tf - solver.stats.time - δt
	new_xf = solver.solver.xf + solver.dxf*δt

	TO.Logging.@info time_remaining
	if solver.opts.live_plotting == :on
		visualize_trajectory_car(solver.solver;save_figure=false)
	end

	reset!(solver.solver, reset_type=:mpc)
    return δt, new_x0, new_xf
end

function need_resolve(solver::DirectGamesSolver)
	# Update g
	cost_expansion(solver.C, solver.obj, solver.Z, solver.model.pu, solver.model.p)

	TO.evaluate!(solver.penalty_constraints, solver.Z)
	TO.jacobian!(solver.penalty_constraints, solver.Z)
	TO.update_active_set!(solver.penalty_constraints, solver.Z)
	penalty_expansion!(solver, solver.Z)

	regularize_primals!(solver.C, solver)

	TO.evaluate!(solver.constraints, solver.Z)
	TO.jacobian!(solver.constraints, solver.Z)
	TO.update_active_set!(solver.constraints, solver.Z)

	TO.evaluate!(solver.dyn_constraints, solver.Z)
	TO.discrete_jacobian!(solver.∇F, solver.model, solver.Z)
	update_g_!(solver)

	# Update cmax
	TO.evaluate!(solver.constraints, solver.Z)
	TO.evaluate!(solver.dyn_constraints, solver.Z)
	TO.max_violation!(solver.constraints)
	TO.max_violation!(solver.dyn_constraints)
	cmax = max(maximum(solver.constraints.c_max),
		maximum(solver.dyn_constraints.c_max))

	if (mean(abs.(solver.g_)) < solver.opts.optimality_constraint_tolerance) && (cmax < solver.opts.constraint_tolerance)
		return false
	end
	return true
end

function generate_new_x0(solver::MPCGamesSolver, δt::T) where T
	# Add noise and propagate dynamics forward
	n,m,N = size(solver)
	new_x0 = mpc_propagate_dynamics(solver.solver.model, solver.solver.Z, δt)

	rnd = SVector{n}(rand(n) .- 0.5)
	new_x0 += solver.opts.noise .* rnd .* δt
	selfish_x0 = Array(solver.solver.x0)
	for (i, ind) in enumerate(solver.opts.selfish_inds)
		selfish_x0[ind] += δt*solver.opts.selfish_dx[i]
	end
	for i in setdiff(1:n, solver.opts.selfish_inds)
		selfish_x0[i] = new_x0[i]
	end
	new_x0 = SVector{n}(selfish_x0)
	return new_x0
end

function mpc_propagate_dynamics(model::L, Z::TO.Traj, Δt::T) where {L,T}
	N = length(Z)
	H = sum([Z[k].dt for k=1:N-1])
	if H < Δt
		@warn "The horizon of the problem is shorter than the queried Δt."
	end
	t = Δt
	k = 0
	x = TO.state(Z[1])
	while t > 0 && k < N
		k +=1
		dt = min(Z[k].dt, t)
		x = discrete_dynamics(
		    TO.DEFAULT_Q,
		    model,
		    TO.state(Z[k]),
		    TO.control(Z[k]),
		    0.0,
		    dt)
		t -= dt
	end

	# k indicates the number of controls that have been applied
	# the last one might not have been used for its full dt.
	return x
end

"""
Stash iteration statistics
"""
function record_iteration!(solver::MPCGamesSolver, δt; resolve::Bool=true)
    solver.stats.iterations += 1
    i = solver.stats.iterations::Int
	solver.stats.time += δt
	solver.stats.solve_time[i] = δt
	x0 = TO.KnotPoint(
		TO.state(solver.solver.Z[1]),
		TO.control(solver.solver.Z[1]),
		δt,
		0.)
	push!(solver.stats.x0, x0)

	if resolve
		solver.stats.solver_iterations[i] = solver.solver.stats.iterations_total
		j = solver.solver.stats.iterations
		solver.stats.cmax[i] = solver.solver.stats.cmax[j]
	else
		solver.stats.solver_iterations[i] = 0
		TO.evaluate!(solver.solver.constraints, solver.solver.Z)
		TO.evaluate!(solver.solver.dyn_constraints, solver.solver.Z)
		TO.max_violation!(solver.solver.constraints)
		TO.max_violation!(solver.solver.dyn_constraints)
		cmax = max(maximum(solver.solver.constraints.c_max),
			maximum(solver.solver.dyn_constraints.c_max))
		solver.stats.cmax[i] = cmax
	end

	cmax_collision, cmax_boundary, cmax_bound =	evaluate_constraint_violation(solver.solver, solver.solver.Z)
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
function evaluate_convergence(solver::MPCGamesSolver)
    # Get current iterations
	i = solver.stats.iterations
	t = solver.stats.time
	# Check constraint violation
	if solver.stats.cmax_collision[i] >= 2e-3
		solver.stats.failure_status = true
		TO.Logging.@info "Collision constraint is violated cmax_collision >= 2e-3."
		return true
	end
	# if solver.stats.cmax_boundary[i] >= 2e-3
	# 	solver.stats.failure_status = true
	# 	TO.Logging.@info "Boundary constraint is violated cmax_boundary >= 2e-3."
	# 	return true
	# end
	if solver.stats.cmax_bound[i] >= 1e-2
		solver.stats.failure_status = true
		TO.Logging.@info "Bound constraint is violated cmax_bound >= 1e-2."
		return true
	end
	# Check total iterations
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

function resample!(solver::MPCGamesSolver{T}) where T
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


function evaluate_constraint_violation(solver::DirectGamesSolver, Z::TO.Traj)
	# # Dynamics
	# TO.evaluate!(solver.dyn_constraints, Z)
	# cmax_dynamics = maximum(solver.dyn_constraints.constraints[1].vals[1])
	# Collision, Boundary, Bound
	cmax_collision = 0.
	cmax_boundary = 0.
	cmax_bound = 0.

	for con in solver.constraints
		evaluate!(con.vals, con.con, Z[1:1], [1])
		if typeof(con.con) <: CollisionConstraint
			cmax_collision = max(cmax_collision, maximum(con.vals[1]))
		elseif typeof(con.con) <: BoundaryConstraint
			cmax_boundary = max(cmax_boundary, maximum(con.vals[1]))
		elseif typeof(con.con) <: BoundConstraint
			cmax_bound = max(cmax_bound, maximum(con.vals[1]))
		end
	end
	# return cmax_dynamics, cmax_collision, cmax_boundary, cmax_bound
	return cmax_collision, cmax_boundary, cmax_bound
end
