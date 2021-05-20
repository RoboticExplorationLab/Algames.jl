################################################################################
# Solver Methods
################################################################################

function newton_solve!(prob::GameProblem{KN,n,m,T,SVd,SVx}) where {KN,n,m,T,SVd,SVx}
	model = prob.model
	core = prob.core
	game_con = prob.game_con
	opts = prob.opts

	# Set initial trajectory
	Random.seed!(opts.seed)
	init_traj!(prob.pdtraj; x0=prob.x0, f=opts.f_init, amplitude=opts.amplitude_init, s=opts.shift)
	init_traj!(prob.pdtraj_trial; x0=prob.x0, f=opts.f_init, amplitude=opts.amplitude_init, s=opts.shift)
	init_traj!(prob.Δpdtraj; x0=prob.x0, f=zeros, amplitude=0.0)

	rollout!(RK3, prob.model, prob.pdtraj.pr)
	rollout!(RK3, prob.model, prob.pdtraj_trial.pr)
	# Set the initial penalties
	prob.pen.ρ = SVector{1,T}([opts.ρ_0])
	prob.pen.ρ_trial = SVector{1,T}([opts.ρ_trial])

	# Reset Statistics and constraints
	reset!(prob.stats)
	opts.dual_reset ? reset!(game_con) : nothing
	# Iterative solve
	out = 0
    for k = 1:opts.outer_iter
		out = k
		# plot_traj!(model, prob.pdtraj.pr)
		# Initialize regularization and failed line search count.
		set!(opts.reg, opts.reg_0)
		LS_count = 0
		opts.inner_print ? display_solver_header() : nothing
		#XXX opts.inner_print ? display_condition_header() : nothing
        for l = 1:opts.inner_iter
			set!(opts.reg, opts.reg_0*l^4)
			LS_count, control_flow = inner_iteration(prob, LS_count, k, l)
			LS_count >= 1 || control_flow == :break ? break : nothing
		end

		# Visualize trajectory
		#XXX live_vis && (opts.outer_iter-k)%1==0 ? visualize!(vis, model, pdtraj.q) : nothing
		# Kick out before updating the penalty and duals
		if k == opts.outer_iter || (
			prob.stats.dyn_vio[end].max < opts.ϵ_dyn &&
			prob.stats.con_vio[end].max < opts.ϵ_con &&
			prob.stats.sta_vio[end].max < opts.ϵ_sta &&
			prob.stats.opt_vio[end].max < opts.ϵ_opt)
			break
		end
		# Dual Ascent
		evaluate!(game_con, prob.pdtraj.pr)
		dual_update!(game_con)
		# Increasing Schedule
		prob.pen.ρ = min.(prob.pen.ρ * opts.ρ_increase, opts.ρ_max)
		penalty_update!(game_con)
    end
	record!(prob.stats, core, model, game_con, prob.pdtraj, out)
    return nothing
end

function inner_iteration(prob::GameProblem, LS_count::Int, k::Int, l::Int)
	core = prob.core
	opts = prob.opts
	# plot_traj!(prob.model, prob.pdtraj.pr)

	# Residual
	residual!(prob, prob.pdtraj)
	opts.regularize ? regularize_residual!(core, opts, prob.pdtraj, prob.pdtraj) : nothing # should do nothing since we regularize around pdtraj
	record!(prob.stats, prob.core, prob.model, prob.game_con, prob.pdtraj, k)
	res_norm = norm(core.res, 1)/length(core.res)

	if prob.stats.opt_vio[end].max < opts.ϵ_opt
		return LS_count, :break
	end
	# Residual Jacobian
	residual_jacobian!(prob)
	regularize_residual_jacobian!(prob)

	Δtraj = - \(lu(core.jac), core.res)
	set_traj!(core, prob.Δpdtraj, Δtraj)

	# Line Search
	α, j = line_search(prob, res_norm)
	failed_ls = j == opts.ls_iter
	failed_ls ? LS_count += 1 : LS_count = 0
	update_traj!(prob.pdtraj, prob.pdtraj, α, prob.Δpdtraj)
	Δ = Δ_step(prob.Δpdtraj, α)
	if Δ < opts.Δ_min
		return LS_count, :break
	end

	opts.inner_print ? display_solver_data(k, l, j, Δ, res_norm, opts.reg) : nothing
	#XXX opts.inner_print ? display_condition_data(res) : nothing
	return LS_count, :continue
end

function line_search(prob::GameProblem, res_norm::T) where {T}
	core = prob.core
	opts = prob.opts

	j = 1
	α = 1.0
	while j < opts.ls_iter
		update_traj!(prob.pdtraj_trial, prob.pdtraj, α, prob.Δpdtraj)
		residual!(prob, prob.pdtraj_trial)
		opts.regularize ? regularize_residual!(core, opts, prob.pdtraj_trial, prob.pdtraj) : nothing# should add regularization term to prevent deviation from pdtraj
		res_norm_trial = norm(core.res, 1)/length(core.res)
		if res_norm_trial <= (1.0-α*opts.β)*res_norm
			LS_count = 0
			break
		else
			α *= opts.α_decrease # step size decrease
			j += 1
		end
	end
	return α, j
end



################################################################################
# Iterative Best Response Solver Methods
################################################################################

function ibr_newton_solve!(prob::GameProblem{KN,n,m,T,SVd,SVx};
		ibr_opts::IBROptions=IBROptions()) where {KN,n,m,T,SVd,SVx}


	# Set initial trajectory
	Random.seed!(opts.seed)
	init_traj!(prob.pdtraj; x0=prob.x0, f=opts.f_init, amplitude=opts.amplitude_init, s=opts.shift)
	init_traj!(prob.pdtraj_trial; x0=prob.x0, f=opts.f_init, amplitude=opts.amplitude_init, s=opts.shift)
	init_traj!(prob.Δpdtraj; x0=prob.x0, f=zeros, amplitude=0.0)

	rollout!(RK3, prob.model, prob.pdtraj.pr)
	rollout!(RK3, prob.model, prob.pdtraj_trial.pr)

	for q = 1:ibr_opts.ibr_iter # need a better stopping criterion, based on convergence
		for i = 1:p
			ibr_newton_solve!(prob, ibr_opts.ordering[i])
		end
	end
	return nothing
end

function ibr_newton_solve!(prob::GameProblem{KN,n,m,T,SVd,SVx}, i::Int) where {KN,n,m,T,SVd,SVx}
	model = prob.model
	core = prob.core
	game_con = prob.game_con
	opts = prob.opts

	# Set the initial penalties
	prob.pen.ρ = SVector{1,T}([opts.ρ_0])
	prob.pen.ρ_trial = SVector{1,T}([opts.ρ_trial])

	# Reset Statistics and constraints
	reset!(prob.stats)
	opts.dual_reset ? reset!(game_con) : nothing
	# Iterative solve
	out = 0
    for k = 1:opts.outer_iter
		out = k
		# plot_traj!(model, prob.pdtraj.pr)
		# Initialize regularization and failed line search count.
		set!(opts.reg, opts.reg_0)
		LS_count = 0
		opts.inner_print ? display_solver_header() : nothing
		#XXX opts.inner_print ? display_condition_header() : nothing
        for l = 1:opts.inner_iter
			set!(opts.reg, opts.reg_0*l^4)
			LS_count, control_flow = ibr_inner_iteration(prob, LS_count, k, l, i)
			LS_count >= 1 || control_flow == :break ? break : nothing
		end

		# Visualize trajectory
		#XXX live_vis && (opts.outer_iter-k)%1==0 ? visualize!(vis, model, pdtraj.q) : nothing
		# Kick out before updating the penalty and duals
		if k == opts.outer_iter || (
			prob.stats.dyn_vio[end].max < opts.ϵ_dyn &&
			prob.stats.con_vio[end].max < opts.ϵ_con &&
			prob.stats.sta_vio[end].max < opts.ϵ_sta &&
			prob.stats.opt_vio[end].max < opts.ϵ_opt)
			break
		end
		# Dual Ascent
		evaluate!(game_con, prob.pdtraj.pr)
		dual_update!(game_con)
		# Increasing Schedule
		prob.pen.ρ = min.(prob.pen.ρ * opts.ρ_increase, opts.ρ_max)
		penalty_update!(game_con)
    end
	record!(prob.stats, core, model, game_con, prob.pdtraj, out)
    return nothing
end

function ibr_inner_iteration(prob::GameProblem, LS_count::Int, k::Int, l::Int, i::Int)
	core = prob.core
	opts = prob.opts
	res_mask = residual_mask(prob, i)
	jac_mask = jacobian_mask(prob, i)
	# plot_traj!(prob.model, prob.pdtraj.pr)

	# Residual
	ibr_residual!(prob, prob.pdtraj, i)
	opts.regularize ? regularize_ibr_residual!(core, opts, prob.pdtraj, prob.pdtraj, i) : nothing # should do nothing since we regularize around pdtraj
	record!(prob.stats, prob.core, prob.model, prob.game_con, prob.pdtraj, k)
	res_norm = norm(core.res[res_mask], 1)/length(core.res[res_mask])

	if prob.stats.opt_vio[end].max < opts.ϵ_opt
		return LS_count, :break
	end
	# Residual Jacobian
	ibr_residual_jacobian!(prob, i)
	regularize_ibr_residual_jacobian!(prob, i)

	Δtraj = zeros(prob.probsize.S)
	Δtraj[res_traj] = - \(lu(core.jac[jac_mask]), core.res[res_mask])
	set_traj!(core, prob.Δpdtraj, Δtraj)

	# Line Search
	α, j = line_search(prob, res_norm)
	failed_ls = j == opts.ls_iter
	failed_ls ? LS_count += 1 : LS_count = 0
	update_traj!(prob.pdtraj, prob.pdtraj, α, prob.Δpdtraj)
	Δ = Δ_step(prob.Δpdtraj, α)
	if Δ < opts.Δ_min
		return LS_count, :break
	end

	opts.inner_print ? display_solver_data(k, l, j, Δ, res_norm, opts.reg) : nothing
	#XXX opts.inner_print ? display_condition_data(res) : nothing
	return LS_count, :continue
end

# ibr_residual! DDDDDDDDDD
# regularize_ibr_residual!DDDDDDDDDDDD
# ibr_residual_jacobian!DDDDDDDDDDDDDDD
# regularize_ibr_residual_jacobian!DDDDDDDDDDDDDDDDD
# res_mask = residual_mask(prob, i)
# jac_mask = jacobian_mask(prob, i)
# record needs to be investigated in
prob.core.res

prob.probsize.S
