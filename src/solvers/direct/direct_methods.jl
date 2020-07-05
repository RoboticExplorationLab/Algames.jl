export
	solve!,
	step!,
	record_iteration!,
	record_inner_iteration!,
	evaluate_convergence,
	evaluate_inner_convergence,
	rollout!,
	regularization_update!,
	regularize_primals!,
	max_violation

	# Generic solve methods
"DirectGames solve method (non-allocating)"
function TO.solve!(solver::DirectGamesSolver{T}) where T<:AbstractFloat
    TO.set_verbosity!(solver.opts)
    TO.clear_cache!(solver.opts)

    solver.stats.iterations = 0
    solver.ρ[1] = 0.0
	solver.dρ[1] = 0.0
	solver.η[1] = 0.0 #####

    n,m,pu,p = size(solver.model)
    n,m,N = size(solver)
    # J = Inf
    # _J = TO.get_J.(solver.obj)

    # Initial rollout
    TO.rollout!(solver)
	for k = 1:N
	    solver.Z̄[k].z = solver.Z[k].z
	end
    for i = 1:p
        cost!(solver.obj[i], solver.Z, fill(Array(pu[i]),N))
    end
    # J_prev = sum.(_J)
	J_prev = indiv_cost(solver)
    for i = 1:solver.opts.iterations
		dt = @elapsed begin
			TO.Logging.@info "Solver iteration = ", i
	        J = step!(solver, J_prev)
	        # check for cost blow up
	        if any(J .> solver.opts.max_cost_value)
	            @warn "Cost exceeded maximum cost"
	            return solver
	        end

			dJ = abs.(J .- J_prev)
			J_prev = copy(J)
		end
 		record_iteration!(solver, J, dJ, dt)
        evaluate_convergence(solver) ? break : nothing #######
		penalty_update!(solver)
    end
    return solver
end

"""
Take one step of DirectGames algorithm (non-allocating)
"""
function step!(solver::DirectGamesSolver, J_prev)
	for i = 1:solver.opts.inner_iterations
		# println("inner iteration = ", i)
		cost_expansion(solver.C, solver.obj, solver.Z, solver.model.pu, solver.model.p)
		###
		###
		TO.evaluate!(solver.penalty_constraints, solver.Z)# solver.Z̄)
		TO.jacobian!(solver.penalty_constraints, solver.Z)# solver.Z̄)
		TO.update_active_set!(solver.penalty_constraints, solver.Z)# solver.Z̄)
		penalty_expansion!(solver, solver.Z)# solver.Z̄)
		###
		###
		regularize_primals!(solver.C, solver)

	    TO.evaluate!(solver.constraints, solver.Z)
	    TO.jacobian!(solver.constraints, solver.Z)
		TO.update_active_set!(solver.constraints, solver.Z)

		TO.evaluate!(solver.dyn_constraints, solver.Z)
		TO.discrete_jacobian!(solver.∇F, solver.model, solver.Z)

		update_g_!(solver)
		if evaluate_inner_convergence(solver) && solver.opts.break_inner_loop
			# @show "break early"
			record_inner_iteration!(solver, NaN)
			break
		end
		update_H_!(solver)
		inner_step!(solver)
	end
	J = indiv_cost(solver)
    return J
end


"""
Stash iteration statistics
"""
function record_iteration!(solver::DirectGamesSolver, J, dJ, dt)
    solver.stats.iterations += 1
    i = solver.stats.iterations::Int
	inner_iter = solver.stats.iterations_inner[i]
	solver.stats.iterations_total += inner_iter

    solver.stats.cost[i] = J
    solver.stats.dJ[i] = dJ

	TO.evaluate!(solver.constraints, solver.Z)
	TO.evaluate!(solver.dyn_constraints, solver.Z)
	TO.max_violation!(solver.constraints)
	TO.max_violation!(solver.dyn_constraints)
	solver.stats.cmax[i] = TO.max_violation(solver)
	solver.stats.runtime += dt

    @logmsg TO.InnerLoop :iter value=i
    @logmsg TO.InnerLoop :cost value=J
    @logmsg TO.InnerLoop :dJ   value=dJ
    if solver.opts.verbose
        print_level(InnerLoop)
    end
    return nothing
end


function record_inner_iteration!(solver::DirectGamesSolver, α)
    i = solver.stats.iterations::Int + 1
	solver.stats.iterations_inner[i] += 1
	j = solver.stats.iterations_inner[i]
	solver.stats.optimality_merit[i][j] = mean(abs.(solver.g_))
	solver.stats.optimality_merit_inf[i][j] = norm(solver.g_, Inf)
	solver.stats.α[i][j] = α
	if solver.opts.record_condition
		solver.stats.H_cond[i][j] = cond(Array(solver.H_))
	else
		solver.stats.H_cond[i][j] = NaN
	end
    return nothing
end


# """
# $(SIGNATURES)
# Check convergence conditions for DirectR
# """
function evaluate_convergence(solver::DirectGamesSolver)
    # Get current iterations
    i = solver.stats.iterations
	j = solver.stats.iterations_inner[i]
	cmax = solver.stats.cmax[i]
	optimality_merit = solver.stats.optimality_merit[i][j]
	optimality_merit_inf = solver.stats.optimality_merit_inf[i][j]

    # Check for cost convergence
    # note the dJ > 0 criteria exists to prevent loop exit when forward pass makes no improvement
    # if all(0.0 .< solver.stats.dJ[i]) && all(solver.stats.dJ[i] .< solver.opts.cost_tolerance) ####
	if (mean(abs.(solver.stats.dJ[i])) < solver.opts.cost_tolerance) && (cmax < solver.opts.constraint_tolerance) &&
		(i>=solver.opts.min_iterations)
		TO.Logging.@info "Outer loop converged: cost_tolerance & constraint_tolerance"
		return true
    end

	if (optimality_merit < solver.opts.optimality_constraint_tolerance) &&
		(optimality_merit_inf < solver.opts.optimality_constraint_tolerance_inf) &&
		(cmax < solver.opts.constraint_tolerance) &&
		(i>=solver.opts.min_iterations)
		TO.Logging.@info "Outer loop converged: optimality_merit & constraint_tolerance"
		return true
	end

    # Check total iterations
    if i >= solver.opts.iterations
		TO.Logging.@info "Outer loop converged: iterations"
		return true
    end

	if solver.stats.runtime > solver.opts.timeout
		TO.Logging.@info "Outer loop converged: timeout"
		return true
	end

    return false
end

function evaluate_inner_convergence(solver::DirectGamesSolver)
	i = solver.stats.iterations + 1
	j = solver.stats.iterations_inner[i]
	# Check optimality
	# if mean(abs.(solver.g_)) < solver.opts.optimality_constraint_tolerance
	# We force the solver to take an inner step for each outer step
	# so that we make progress even if the optimality constraint is respected.
	# Indeed, the optimality constraint ~scales with the penalty term
	# so it is small at the beginning.
	if (mean(abs.(solver.g_)) < solver.opts.optimality_constraint_tolerance) &&
		(j>=solver.opts.min_steps_per_iteration)
		@debug "Inner loop converged optimality"
		return true
	end

	# Check total iterations
    if j >= solver.opts.inner_iterations
		@debug "Inner loop converged iterations"
        return true
    end

	# Outer loop update if line search is repeatedly unsuccessful
    if solver.stats.dJ_zero_counter > solver.opts.dJ_counter_limit
		@debug "Inner loop converged dJ zero counter"
		solver.stats.dJ_zero_counter = 0 # reset to 0 before the new outer loop
        return true
    end
    return false
end



"Simulate the forward the dynamics open-loop"
@inline rollout!(solver::DirectGamesSolver) = rollout!(solver.model, solver.Z, solver.x0)

function rollout!(model::AbstractGameModel, Z::Traj, x0)
    Z[1].z = [x0; control(Z[1])]
    for k = 2:length(Z)
        TO.propagate_dynamics(TO.DEFAULT_Q, model, Z[k], Z[k-1])
    end
end


# """
# $(SIGNATURES)
# Update the regularzation for the DirectR backward pass
# """
function regularization_update!(solver::DirectGamesSolver,status::Symbol=:increase)
    if status == :increase # increase regularization
        # @logmsg InnerLoop "Regularization Increased"
        solver.dρ[1] = max(solver.dρ[1]*solver.opts.bp_reg_increase_factor, solver.opts.bp_reg_increase_factor)
		solver.ρ[1] = max(solver.ρ[1]*solver.dρ[1], solver.opts.bp_reg_min)
		# solver.η[1] = max(solver.η[1]*solver.dρ[1], solver.opts.bp_reg_min)
        # if solver.ρ[1] > solver.opts.bp_reg_max
        #     @warn "Max regularization exceeded"
        # end
    elseif status == :decrease # decrease regularization
        # TODO: Avoid divides by storing the decrease factor (divides are 10x slower)
        solver.dρ[1] = min(solver.dρ[1]/solver.opts.bp_reg_increase_factor, 1.0/solver.opts.bp_reg_increase_factor)
		solver.ρ[1] = solver.ρ[1]*solver.dρ[1]*(solver.ρ[1]*solver.dρ[1]>solver.opts.bp_reg_min)
		# solver.η[1] = solver.η[1]*solver.dρ[1]*(solver.η[1]*solver.dρ[1]>solver.opts.bp_reg_min)
    end
end


function regularize_primals!(C::Vector{E}, solver::TO.AbstractSolver) where {T, E<:TO.CostExpansion}
	n,m,N = size(solver)
	n,m,pu,p = size(solver.model)
	η = solver.η
	ηx = η[1]*Diagonal(@SVector ones(n))
	for i = 1:p
		ηu = η[1]*Diagonal(@SVector ones(length(pu[i])))
		for k = 1:N-1
			C[i].xx[k] += ηx
			C[i].uu[k] += ηu
		end
		C[i].xx[N] += ηx
	end
	return nothing
end


function TO.max_violation(solver::DirectGamesSolver)
	cmax = 0.
	if !isempty(solver.constraints.constraints)
		cmax = max(maximum(solver.constraints.c_max),
			maximum(solver.dyn_constraints.c_max))
	end
	return cmax
end
