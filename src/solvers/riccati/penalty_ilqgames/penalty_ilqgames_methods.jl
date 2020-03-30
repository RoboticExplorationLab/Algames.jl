export
	solve!,
	set_penalty!,
	step!,
	penalty_expansion!,
	record_iteration!,
	gradient_todorov!,
	evaluate_convergence,
	nash_feedback_backwardpass!,
	nash_open_loop_backwardpass!,
	forwardpass!,
	rollout!,
	regularization_update!,
	max_violation

# Generic solve methods
"PenaltyiLQGames solve method (non-allocating)"
function TO.solve!(solver::PenaltyiLQGamesSolver{T}) where T<:AbstractFloat
    TO.set_verbosity!(solver.opts)
    TO.clear_cache!(solver.opts)

    solver.stats.iterations = 0
    solver.ρ[1] = 0.0
    solver.dρ[1] = 0.0
    # reset!(solver)
    Z = solver.Z; Z̄ = solver.Z̄;

    n,m,pu,p = size(solver.model)
    n,m,N = size(solver)
    J = Inf
    _J = TO.get_J.(solver.obj)


    # Initial rollout
    rollout!(solver)

    for i = 1:p
        cost!(solver.obj[i], solver.Z, fill(pu[i],N))
    end
    for i = 1:p
        for con in solver.constraints.constraints
            TO.cost!(_J[i], con, solver.Z)
        end
    end
    J_prev = sum.(_J)


    set_penalty!(solver.constraints, solver.pen)
    for i = 1:solver.opts.iterations
        TO.Logging.@info "Solver iteration = ", i
        J, ΔV, α = step!(solver, J_prev)

		# check for cost blow up
        if any(J .> solver.opts.max_cost_value)
            @warn "Cost exceeded maximum cost"
            return solver
        end

        for k = 1:N
            Z[k].z = Z̄[k].z
        end

        dJ = abs.(J .- J_prev)
        J_prev = copy(J)
        gradient_todorov!(solver)

        record_iteration!(solver, J, dJ, ΔV, α)
        evaluate_convergence(solver) ? break : nothing
    end
    return solver
end


"""
Take one step of PenaltyiLQGames algorithm (non-allocating)
"""
function step!(solver::PenaltyiLQGamesSolver, J)
    # Constraints
    TO.evaluate!(solver.constraints, solver.Z)
    TO.jacobian!(solver.constraints, solver.Z)
    TO.update_active_set!(solver.constraints, solver.Z)

    TO.discrete_jacobian!(solver.∇F, solver.model, solver.Z)
    cost_expansion(solver.C, solver.obj, solver.Z, solver.model.pu, solver.model.p)
    penalty_expansion!(solver, solver.Z)
    if solver.opts.eq_type == :nash && solver.opts.info_pattern == :feedback
        ΔV = nash_feedback_backwardpass!(solver)
        J, α = forwardpass!(solver, ΔV, J)
        return J, ΔV, α
    elseif solver.opts.eq_type == :nash && solver.opts.info_pattern == :open_loop
        ΔV = nash_open_loop_backwardpass!(solver)
        J, α = forwardpass!(solver, ΔV, J)
        return J, ΔV, α
    else
        @warn "This type of equilibrium is not supported yet."
    end
end

function penalty_expansion!(solver::PenaltyiLQGamesSolver{T}, Z::TO.Traj) where T
    C = solver.C
    conSet = solver.constraints
    n,m,pu,p = size(solver.model)
    n,m,N = size(solver)
    for i = 1:p
        for con in conSet.constraints
            if typeof(con).parameters[2] == State
                for (j,k) in enumerate(intersect(con.inds,1:N))
                    Iμ = TO.penalty_matrix(con,j)
                    C[i].xx[k] += con.∇c[j]'*Iμ*con.∇c[j]
                    C[i].x[k] += con.∇c[j]'*Iμ*con.vals[j]
                end
            elseif typeof(con).parameters[2] == Control
                for (j,k) in enumerate(intersect(con.inds,1:N-1))
                    Iμ = TO.penalty_matrix(con,j)
                    C[i].uu[k] += con.∇c[j]'*Iμ*con.∇c[j]
                    C[i].u[k] += con.∇c[j]'*Iμ*con.vals[j]
                end
			else
				@warn(typeof(con), " this type of constraint is not supported.")
            end
        end
    end
    return nothing
end


"""
Stash iteration statistics
"""
function record_iteration!(solver::PenaltyiLQGamesSolver, J, dJ, ΔV, α)
    solver.stats.iterations += 1
    i = solver.stats.iterations::Int
    solver.stats.cost[i] = J
    solver.stats.dJ[i] = dJ
    solver.stats.gradient[i] = mean(solver.grad)

    TO.evaluate!(solver.constraints, solver.Z)
	TO.max_violation!(solver.constraints)
	solver.stats.cmax[i] = TO.max_violation(solver)
	solver.stats.ΔV[i] = Array.(ΔV)
	solver.stats.α[i] = α

    @logmsg TO.InnerLoop :iter value=i
    @logmsg TO.InnerLoop :cost value=J
    @logmsg TO.InnerLoop :dJ   value=dJ
    @logmsg TO.InnerLoop :grad value=solver.stats.gradient[i]
    if solver.opts.verbose
        print_level(InnerLoop)
    end
    return nothing
end

# """
# $(SIGNATURES)
#     Calculate the problem gradient using heuristic from iLQG (Todorov) solver
# """
function gradient_todorov!(solver::PenaltyiLQGamesSolver)
    for k in eachindex(solver.d)
		solver.grad[k] = mean( abs.(solver.d[k]) ./ (abs.(control(solver.Z[k])) .+ 1) ) ################################
		# solver.grad[k] = maximum( abs.(solver.d[k]) ./ (abs.(control(solver.Z[k])) .+ 1) )
    end
end


# """
# $(SIGNATURES)
# Check convergence conditions for iLQR
# """
function evaluate_convergence(solver::PenaltyiLQGamesSolver)
    # Get current iterations
    i = solver.stats.iterations

    # Check for cost convergence
    # note the dJ > 0 criteria exists to prevent loop exit when forward pass makes no improvement
	if mean(abs.(solver.stats.dJ[i])) < solver.opts.cost_tolerance && (solver.stats.cmax[i] < solver.opts.constraint_tolerance)
        TO.Logging.@info "dJ tolerance < cost_tolerance"
        return true
    end

    # Check for gradient convergence
    if solver.stats.gradient[i] < solver.opts.gradient_norm_tolerance && (solver.stats.cmax[i] < solver.opts.constraint_tolerance)
        TO.Logging.@info "gradient tolerance < gradient_norm_tolerance"
        return true
    end

    # Check total iterations
    if i >= solver.opts.iterations
        TO.Logging.@info "iterations convergence"
        return true
    end

    return false
end

# """
# $(SIGNATURES)
# Calculates the optimal feedback gains K,d as well as the 2nd Order approximation of the
# Cost-to-Go, using a backward Riccati-style recursion. (non-allocating)
# """
function nash_feedback_backwardpass!(solver::PenaltyiLQGamesSolver{T,QUAD}) where {T,QUAD<:QuadratureRule}
    n,m,N = size(solver)
    n,m,pu,p = size(solver.model)

    # Objective
    obj = solver.obj
    model = solver.model


    # Extract variables
    Z = solver.Z; K = solver.K; d = solver.d;
    S = solver.S
    Q = solver.Q
    C = solver.C

    # Terminal cost-to-go
    for i = 1:p
        S[i].xx[N] = C[i].xx[N]
        S[i].x[N] = C[i].x[N]
    end

    # Initialize expected change in cost-to-go
    ΔV = [@SVector zeros(2) for i=1:p]

    k = N-1
    while k > 0
        ix = Z[k]._x
        iu = Z[k]._u
        x = state(Z[k])
        u = control(Z[k])

        fdx = solver.∇F[k][ix,ix] # Ak [n,n]
        fdu = solver.∇F[k][ix,iu] # Bk [n,m]

        for i = 1:p
            Q[i].x[k] = fdx'*S[i].x[k+1]
            Q[i].u[k] = fdu'*S[i].x[k+1]
            Q[i].xx[k] = fdx'S[i].xx[k+1]*fdx
            Q[i].uu[k] = fdu'*S[i].xx[k+1]*fdu
            Q[i].ux[k] = fdu'*S[i].xx[k+1]*fdx
        end
        for i = 1:p
            H = Hermitian(Array((C[i].uu[k] + Q[i].uu[k][pu[i], pu[i]])))
            if !isposdef(H)
                @warn "Player $i's value function 2nd order approximation is not positive definite"
            end
        end

        # Regularization
        for i = 1:p
            if solver.opts.bp_reg_type == :state
                Quu_reg = Q[i].uu[k] + solver.ρ[1]*fdu'fdu
                Qux_reg = Q[i].ux[k] + solver.ρ[1]*fdu'fdx
            elseif solver.opts.bp_reg_type == :control
                Quu_reg = Q[i].uu[k] + solver.ρ[1]*I
                Qux_reg = Q[i].ux[k]
            end
            solver.N_[pu[i],:] = Quu_reg[:,pu[i]]'
            solver.N_[pu[i],pu[i]] += C[i].uu[k]' # gik_uiui^T
            solver.M_[pu[i],:] = - Qux_reg[pu[i],:]
            solver.M_[pu[i],:] += - C[i].ux[k] # gik_uix

            solver.m_[pu[i]] = - Q[i].u[k][pu[i]]
            solver.m_[pu[i]] += - C[i].u[k] # gik_ui
        end

        # Compute gains
        K[k] = solver.N_\solver.M_ # [m,n]
        d[k] = solver.N_\solver.m_ # [m]
        # Calculate cost-to-go (using unregularized Cuu Quu and Cux Qux)
        for i = 1:p
            S[i].xx[k] = Q[i].xx[k] +
                K[k]'*Q[i].uu[k]*K[k] +
                Q[i].ux[k]'*K[k] +
                K[k]'*Q[i].ux[k]
            S[i].xx[k] += C[i].xx[k] + # gik_xx
                K[k][pu[i],:]' * C[i].uu[k] * K[k][pu[i],:] + # gik_uiui
                C[i].ux[k]' * K[k][pu[i],:] + # gik_uix
                K[k][pu[i],:]' * C[i].ux[k] # gik_uix
            S[i].xx[k] = 0.5*(S[i].xx[k] + S[i].xx[k]')

            S[i].x[k] = Q[i].x[k] +
                K[k]'* Q[i].u[k] +
                Q[i].ux[k]'* d[k] +
                K[k]'* Q[i].uu[k] * d[k]
            S[i].x[k] += C[i].x[k] + # gik_x
                K[k][pu[i],:]'* C[i].u[k] + # gik_ui
                C[i].ux[k]'* d[k][pu[i]] + # gik_uix
                K[k][pu[i],:]'* C[i].uu[k]'* d[k][pu[i]] # gik_uiui
        end

        # calculated change is cost-to-go over entire trajectory
        ΔV += [@SVector [d[k]'*Q[i].u[k], 0.5*d[k]'*Q[i].uu[k]*d[k]] for i=1:p]
        k -= 1
    end
    regularization_update!(solver, :decrease)
    return ΔV

end


function nash_open_loop_backwardpass!(solver::PenaltyiLQGamesSolver{T,QUAD}) where {T,QUAD<:QuadratureRule}
    n,m,N = size(solver)
    n,m,pu,p = size(solver.model)

    # Objective
    obj = solver.obj
    model = solver.model

    # Extract variables
    Z = prob.Z; K = solver.K; d = solver.d;
    S = solver.S
    Q = solver.Q
    C = solver.C

    # Terminal cost-to-go
    for i = 1:p
        S[i].xx[N] = C[i].xx[N] # αiN
        S[i].x[N] = C[i].x[N] # βiN
    end

    # Initialize expected change in cost-to-go
    ΔV = [@SVector zeros(2) for i=1:p]

    k = N-1
    while k > 0
        # println("backwardpass step = ", k)
        ix = Z[k]._x
        iu = Z[k]._u
        x = state(Z[k])
        u = control(Z[k])

        fdx = solver.∇F[k][ix,ix] # Ak [n,n]
        fdu = solver.∇F[k][ix,iu] # Bk [n,m]

        # Regularization
        for i = 1:p
            ##### if solver.opts.bp_reg_type == :state
            #     Quu_reg = Q[i].uu[k] + solver.ρ[1]*fdu'fdu
            #     Qux_reg = Q[i].ux[k] + solver.ρ[1]*fdu'fdx
            # elseif solver.opts.bp_reg_type == :control
            #     Quu_reg = Q[i].uu[k] + solver.ρ[1]*I
            #     Qux_reg = Q[i].ux[k]
            # end
            solver.N_[pu[i],:] = fdu[:,pu[i]]'*S[i].xx[k+1]*fdu
            solver.N_[pu[i],pu[i]] += C[i].uu[k] # gik_uiui
            solver.M_[pu[i],:] = - fdu[:,pu[i]]'*S[i].xx[k+1]*fdx
            solver.M_[pu[i],:] += - C[i].ux[k] # gik_uix

            solver.m_[pu[i]] = - fdu[:,pu[i]]'*S[i].x[k+1]
            solver.m_[pu[i]] += - C[i].u[k] # gik_ui
        end

        for i = 1:p
            if !isposdef(Array(solver.N_[pu[i],pu[i]]))
                @warn "Player $i 's value function 2nd order approximation is not positive definite"
            end
        end

        # Compute gains
        K[k] = solver.N_\solver.M_ # [m,n]
        d[k] = solver.N_\solver.m_ # [m]

        for i = 1:p
            Q[i].xx[k] = S[i].xx[k+1]*(fdx+fdu*K[k]) # Pik
            Q[i].x[k] = S[i].xx[k+1]*fdu*d[k] + S[i].x[k+1] # pik
        end

        # Calculate cost-to-go (using unregularized Cuu Quu and Cux Qux)
        for i = 1:p
            S[i].xx[k] = C[i].xx[k] + # gik_xx
                C[i].ux[k]'* K[k][pu[i],:] + # gik_xui
                fdx'* Q[i].xx[k]
            S[i].x[k] = C[i].ux[k]'* d[k][pu[i]] + # gik_xui
                fdx'* Q[i].x[k] + # gik_ui
                C[i].x[k] # gik_x
        end

        #### # calculated change is cost-to-go over entire trajectory
        #### ΔV += [@SVector [d[k]'*Q[i].u[k], 0.5*d[k]'*Q[i].uu[k]*d[k]] for i=1:p]
        k -= 1
    end
    regularization_update!(solver, :decrease)
    return ΔV

end



# """
# $(SIGNATURES)
# Simulate the system forward using the optimal feedback gains from the backward pass,
# projecting the system on the dynamically feasible subspace. Performs a line search to ensure
# adequate progress on the nonlinear problem.
# """
function forwardpass!(solver::PenaltyiLQGamesSolver{T}, ΔV, J_prev) where T
    # println("ΔV = ", [round.(ΔVi, digits=3) for ΔVi in ΔV])
	n,m,pu,p = size(solver.model)
	n,m,N = size(solver)
    Z = solver.Z; Z̄ = solver.Z̄
    obj = solver.obj

    _J = TO.get_J.(obj)
    J::Vector{Float64} = [Inf for i=1:p]
    α = 1.0
    iter = 0
    z = -1000.0
    expected = 0.0
    flag = true
    one_pass = true
    while (z ≤ solver.opts.line_search_lower_bound || z > solver.opts.line_search_upper_bound) #### && J >= J_prev
        # Check that maximum number of line search decrements has not occured
        if iter > solver.opts.iterations_linesearch
            TO.Logging.@info "iter > solver.opts.iterations_linesearch"
            for k in eachindex(Z)
                Z̄[k].z = Z[k].z
            end
            for i = 1:p
                cost!(obj[i], Z̄, fill(pu[i],N))
            end
            # cost!(obj, Z̄)
            for i = 1:p
                for con in solver.constraints.constraints
                    TO.cost!(_J[i], con, solver.Z̄)
                end
            end
            J = sum.(_J)

            z = 0
            α = 0.0
            # flag = rollout!(solver, α)#####################
			α /= 2.0
            # @warn "we changed the α"
            expected = 0.0

            regularization_update!(solver, :increase)
            solver.ρ[1] += solver.opts.bp_reg_fp
            break
        end

        # Otherwise, rollout a new trajectory for current alpha
        flag = rollout!(solver, α)

        # Check if rollout completed
        if ~flag
            TO.Logging.@info "flagged"
            # Reduce step size if rollout returns non-finite values (NaN or Inf)
            # @logmsg InnerIters "Non-finite values in rollout"
            iter += 1
            α /= 2.0
            continue
        end

        # Calcuate cost
        for i = 1:p
            cost!(obj[i], Z̄, fill(pu[i],N))
        end
        # cost!(obj, Z̄)
        for i = 1:p
            for con in solver.constraints.constraints
                TO.cost!(_J[i], con, solver.Z̄)
            end
        end
        J = sum.(_J)

		norm_Z = zeros(N)
		norm_Δ = zeros(N)
		for k = 1:N
			norm_Z[k] = norm(solver.Z[k].z, 1)
			norm_Δ[k] = norm(solver.Z[k].z - solver.Z̄[k].z, 1)
		end
        expected = [-α*(ΔV[i][1] + α*ΔV[i][2]) for i=1:p]
		# z::Float64 = mean((J_prev .- J)./expected) ####\
		# z::Float64 = mean((J_prev .- J)./J_prev) ####\
		z::T = sum(norm_Δ) / sum(norm_Z)
		# @show z
        iter += 1
        α /= 2.0
    end

    @logmsg TO.InnerLoop :expected value=expected
    @logmsg TO.InnerLoop :z value=z
    @logmsg TO.InnerLoop :α value=2*α
    @logmsg TO.InnerLoop :ρ value=solver.ρ[1]

    return J, 2*α
end



# """
# $(SIGNATURES)
# Simulate forward the system with the optimal feedback gains from the iLQR backward pass.
# (non-allocating)
# """
function rollout!(solver::PenaltyiLQGamesSolver{T,Q}, α) where {T,Q}
    Z = solver.Z; Z̄ = solver.Z̄
    K = solver.K; d = solver.d;

    Z̄[1].z = [solver.x0; control(Z[1])]

    temp = 0.0

    for k = 1:solver.N-1
        δx = TO.state_diff(solver.model, state(Z̄[k]), state(Z[k]))
        ū = control(Z[k]) + K[k]*δx + α*d[k]
        TO.set_control!(Z̄[k], ū)

        Z̄[k+1].z = [TO.discrete_dynamics(Q, solver.model, Z̄[k]);
            control(Z[k+1])]


        temp = norm(Z̄[k+1].z)
        if temp > solver.opts.max_state_value
            return false
        end
    end
    return true
end


"Simulate the forward the dynamics open-loop"
@inline rollout!(solver::PenaltyiLQGamesSolver) = TO.rollout!(solver.model, solver.Z, solver.x0)

# """
# $(SIGNATURES)
# Update the regularzation for the iLQR backward pass
# """
function regularization_update!(solver::PenaltyiLQGamesSolver,status::Symbol=:increase)
    if status == :increase # increase regularization
        # @logmsg InnerLoop "Regularization Increased"
        solver.dρ[1] = max(solver.dρ[1]*solver.opts.bp_reg_increase_factor, solver.opts.bp_reg_increase_factor)
        solver.ρ[1] = max(solver.ρ[1]*solver.dρ[1], solver.opts.bp_reg_min)
        # if solver.ρ[1] > solver.opts.bp_reg_max
        #     @warn "Max regularization exceeded"
        # end
    elseif status == :decrease # decrease regularization
        # TODO: Avoid divides by storing the decrease factor (divides are 10x slower)
        solver.dρ[1] = min(solver.dρ[1]/solver.opts.bp_reg_increase_factor, 1.0/solver.opts.bp_reg_increase_factor)
        solver.ρ[1] = solver.ρ[1]*solver.dρ[1]*(solver.ρ[1]*solver.dρ[1]>solver.opts.bp_reg_min)
    end
end


function TO.max_violation(solver::PenaltyiLQGamesSolver)
	cmax = 0.
	if !isempty(solver.constraints.constraints)
		cmax = 	maximum(solver.constraints.c_max)
	end
	return cmax
end
