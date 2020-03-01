export
    solve!,
    step!,
    record_iteration!,
    gradient_todorov!,
    evaluate_convergence,
    nash_feedback_backwardpass!,
    nash_open_loop_backwardpass!,
    forwardpass!,
    rollout!,
    regularization_update!

# Generic solve methods
"iLQGames solve method (non-allocating)"
function TO.solve!(solver::iLQGamesSolver{T}) where T<:AbstractFloat
    TO.set_verbosity!(solver.opts)
    TO.clear_cache!(solver.opts)

    solver.stats.iterations = 0
    solver.ρ[1] = 0.0
    solver.dρ[1] = 0.0
    # reset!(solver)
    # to = solver.stats[:timer]
    Z = solver.Z; Z̄ = solver.Z̄;

    n,m,pl,p = size(solver.model)
    n,m,N = size(solver)
    J = Inf
    _J = TO.get_J.(solver.obj)

    # Initial rollout
    rollout!(solver)

    for i = 1:p
        # cost!(solver.obj[i], prob.Z, fill(pl[i],N)) #### works
        cost!(solver.obj[i], solver.Z, fill(pl[i],N))
    end
    # cost!(solver.obj, Z)
    J_prev = sum.(_J)

    for i = 1:solver.opts.iterations
        println("Solver iteration = ", i)
        J = step!(solver, J_prev)

        # check for cost blow up
        if any(J .> solver.opts.max_cost_value)
            # @warn "Cost exceeded maximum cost"
            return solver
        end

        for k = 1:N
            Z[k].z = Z̄[k].z
        end

        dJ = abs.(J .- J_prev)
        J_prev = copy(J)
        gradient_todorov!(solver)

        record_iteration!(solver, J, dJ)
        evaluate_convergence(solver) ? break : nothing
    end
    return solver
end


"""
Take one step of iLQGames algorithm (non-allocating)
"""
function step!(solver::iLQGamesSolver, J)
    Z = solver.Z
    # TO.state_diff_jacobian!(solver.G, solver.model, Z)
    TO.discrete_jacobian!(solver.∇F, solver.model, Z)
    cost_expansion(solver.C, solver.obj, solver.Z, solver.model.pl, solver.model.p)
    # println("before barckward")
    if solver.opts.eq_type == :nash && solver.opts.info_pattern == :feedback
        # println("before feedback barckward")
        ΔV = nash_feedback_backwardpass!(solver)
        J = forwardpass!(solver, ΔV, J)
        # println("J after feedback FP = ", J)
        return J
    elseif solver.opts.eq_type == :nash && solver.opts.info_pattern == :open_loop
        # println("before openloop barckward")
        ΔV = nash_open_loop_backwardpass!(solver)
        J = forwardpass!(solver, ΔV, J)
        # println("J after openloop FP = ", J)
        return J
    else
        @warn "This type of equilibrium is not supported yet."
    end
end

"""
Stash iteration statistics
"""
function record_iteration!(solver::iLQGamesSolver, J, dJ)
    solver.stats.iterations += 1
    i = solver.stats.iterations::Int
    solver.stats.cost[i] = J
    solver.stats.dJ[i] = dJ
    solver.stats.gradient[i] = mean(solver.grad)

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
function gradient_todorov!(solver::iLQGamesSolver)
    for k in eachindex(solver.d)
        solver.grad[k] = maximum( abs.(solver.d[k]) ./ (abs.(control(solver.Z[k])) .+ 1) )
    end
end


# """
# $(SIGNATURES)
# Check convergence conditions for iLQR
# """
function evaluate_convergence(solver::iLQGamesSolver)
    # Get current iterations
    i = solver.stats.iterations

    # Check for cost convergence
    # note the dJ > 0 criteria exists to prevent loop exit when forward pass makes no improvement
    # if all(0.0 .< solver.stats.dJ[i]) && all(solver.stats.dJ[i] .< solver.opts.cost_tolerance) ####
    if mean(abs.(solver.stats.dJ[i])) < solver.opts.cost_tolerance
        @show solver.stats.dJ[i]
        println("dJ tolerance < cost_tolerance")
        return true
    end

    # Check for gradient convergence
    if solver.stats.gradient[i] < solver.opts.gradient_norm_tolerance
        println("gradient tolerance < gradient_norm_tolerance")
        return true
    end

    # Check total iterations
    if i >= solver.opts.iterations
        println("iterations convergence")
        return true
    end

    # Outer loop update if forward pass is repeatedly unsuccessful
    if solver.stats.dJ_zero_counter > solver.opts.dJ_counter_limit ####
        println("dJ_zero_counter convergence")
        return true
    end

    return false
end

# """
# $(SIGNATURES)
# Calculates the optimal feedback gains K,d as well as the 2nd Order approximation of the
# Cost-to-Go, using a backward Riccati-style recursion. (non-allocating)
# """
function nash_feedback_backwardpass!(solver::iLQGamesSolver{T,QUAD}) where {T,QUAD<:QuadratureRule}
    # println("****************we are in N fB BP!")
    n,m,N = size(solver)
    n,m,pl,p = size(solver.model)

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
        # println("backwardpass step = ", k)
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
            if !isposdef(Array(C[i].uu[k] + Q[i].uu[k][pl[i], pl[i]]))
                @warn "Player $i 's value function 2nd order approximation is not positive definite"
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
            solver.N_[pl[i],:] = Quu_reg[:,pl[i]]'
            solver.N_[pl[i],pl[i]] += C[i].uu[k]' # gik_uiui^T
            solver.M_[pl[i],:] = - Qux_reg[pl[i],:]
            solver.M_[pl[i],:] += - C[i].ux[k] # gik_uix

            solver.m_[pl[i]] = - Q[i].u[k][pl[i]]
            solver.m_[pl[i]] += - C[i].u[k] # gik_ui
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
                K[k][pl[i],:]' * C[i].uu[k] * K[k][pl[i],:] + # gik_uiui
                C[i].ux[k]' * K[k][pl[i],:] + # gik_uix
                K[k][pl[i],:]' * C[i].ux[k] # gik_uix
            S[i].xx[k] = 0.5*(S[i].xx[k] + S[i].xx[k]')

            S[i].x[k] = Q[i].x[k] +
                K[k]'* Q[i].u[k] +
                Q[i].ux[k]'* d[k] +
                K[k]'* Q[i].uu[k] * d[k]
            S[i].x[k] += C[i].x[k] + # gik_x
                K[k][pl[i],:]'* C[i].u[k] + # gik_ui
                C[i].ux[k]'* d[k][pl[i]] + # gik_uix
                K[k][pl[i],:]'* C[i].uu[k]'* d[k][pl[i]] # gik_uiui
        end

        # calculated change is cost-to-go over entire trajectory
        ΔV += [@SVector [d[k]'*Q[i].u[k], 0.5*d[k]'*Q[i].uu[k]*d[k]] for i=1:p]
        # println("+= ΔV = ", [@SVector [d[k]'*Q[i].u[k], 0.5*d[k]'*Q[i].uu[k]*d[k]] for i=1:p])
        k -= 1
    end
    # println("end of backward pass ΔV = ", ΔV)
    regularization_update!(solver, :decrease)
    return ΔV

end


function nash_open_loop_backwardpass!(solver::iLQGamesSolver{T,QUAD}) where {T,QUAD<:QuadratureRule}
    n,m,N = size(solver)
    n,m,pl,p = size(solver.model)

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
            solver.N_[pl[i],:] = fdu[:,pl[i]]'*S[i].xx[k+1]*fdu
            solver.N_[pl[i],pl[i]] += C[i].uu[k] # gik_uiui
            solver.M_[pl[i],:] = - fdu[:,pl[i]]'*S[i].xx[k+1]*fdx
            solver.M_[pl[i],:] += - C[i].ux[k] # gik_uix

            solver.m_[pl[i]] = - fdu[:,pl[i]]'*S[i].x[k+1]
            solver.m_[pl[i]] += - C[i].u[k] # gik_ui
        end

        for i = 1:p
            if !isposdef(Array(solver.N_[pl[i],pl[i]]))
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
                C[i].ux[k]'* K[k][pl[i],:] + # gik_xui
                fdx'* Q[i].xx[k]
            S[i].x[k] = C[i].ux[k]'* d[k][pl[i]] + # gik_xui
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
function forwardpass!(solver::iLQGamesSolver, ΔV, J_prev)
    println("ΔV = ", [round.(ΔVi, digits=3) for ΔVi in ΔV])

    Z = solver.Z; Z̄ = solver.Z̄
    obj = solver.obj

    _J = TO.get_J.(obj)
    J::Vector{Float64} = fill(Inf,p)
    α = 1.0
    iter = 0
    z = -1.0
    expected = 0.0
    flag = true
    one_pass = true
    while (z ≤ solver.opts.line_search_lower_bound || z > solver.opts.line_search_upper_bound) #### && J >= J_prev
    # while one_pass#### && J >= J_prev
    #     one_pass = false
        # println("before z = ", z)
        # Check that maximum number of line search decrements has not occured
        if iter > solver.opts.iterations_linesearch
            println("iter > solver.opts.iterations_linesearch")
            for k in eachindex(Z)
                Z̄[k].z = Z[k].z
            end
            for i = 1:p
                cost!(obj[i], Z̄, fill(pl[i],N))
            end
            # cost!(obj, Z̄)
            J = sum.(_J)

            z = 0
            α = 0.0
            expected = 0.0

            regularization_update!(solver, :increase)
            solver.ρ[1] += solver.opts.bp_reg_fp
            break
        end

        # Otherwise, rollout a new trajectory for current alpha
        flag = rollout!(solver, α)

        # Check if rollout completed
        if ~flag
            println("flagged")
            # Reduce step size if rollout returns non-finite values (NaN or Inf)
            # @logmsg InnerIters "Non-finite values in rollout"
            iter += 1
            α /= 2.0
            continue
        end

        # Calcuate cost
        for i = 1:p
            cost!(obj[i], Z̄, fill(pl[i],N))
        end
        # cost!(obj, Z̄)
        J = sum.(_J)
        expected = [-α*(ΔV[i][1] + α*ΔV[i][2]) for i=1:p]
        # @show mean((J_prev .- J)./expected)
        z::Float64 = mean((J_prev .- J)./expected) ####
        # if expected > 0.0
        #     zs = (J_prev .- J)./expected
        #     z::Float64 = mean(zs) #### need to understand this rule and adpat it to the non cooperative case
        # else
        #     z = -1.0
        # end
        # println("after  z = ", z)

        iter += 1
        α /= 2.0
    end

    # if J > J_prev
    #     error("Error: Cost increased during Forward Pass")
    # end
    println("α = ", α*2.0)

    @logmsg TO.InnerLoop :expected value=expected
    @logmsg TO.InnerLoop :z value=z
    @logmsg TO.InnerLoop :α value=2*α
    @logmsg TO.InnerLoop :ρ value=solver.ρ[1]

    return J
end



# """
# $(SIGNATURES)
# Simulate forward the system with the optimal feedback gains from the iLQR backward pass.
# (non-allocating)
# """
function rollout!(solver::iLQGamesSolver{T,Q}, α) where {T,Q}
    Z = solver.Z; Z̄ = solver.Z̄
    K = solver.K; d = solver.d;

    Z̄[1].z = [solver.x0; control(Z[1])]

    temp = 0.0

    for k = 1:solver.N-1
        δx = TO.state_diff(solver.model, state(Z̄[k]), state(Z[k]))
        ū = control(Z[k]) + K[k]*δx + α*d[k]
        TO.set_control!(Z̄[k], ū)

        # Z̄[k].z = [state(Z̄[k]); control(Z[k]) + δu]
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
@inline rollout!(solver::iLQGamesSolver) = rollout!(solver.model, solver.Z, solver.x0)

# """
# $(SIGNATURES)
# Update the regularzation for the iLQR backward pass
# """
function regularization_update!(solver::iLQGamesSolver,status::Symbol=:increase)
    # println("reg $(status)")
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
