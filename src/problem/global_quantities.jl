################################################################################
# Residual
################################################################################
function residual!(prob::GameProblem{KN,n,m,T,SVd,SVx}) where {KN,n,m,T,SVd,SVx}
	residual!(prob, prob.pdtraj)
	return nothing
end

function residual!(prob::GameProblem{KN,n,m,T,SVd,SVx}, pdtraj::PrimalDualTraj{KN,n,m,T,SVd}) where {KN,n,m,T,SVd,SVx}
	N = prob.probsize.N
	p = prob.probsize.p
    pu = prob.probsize.pu
	core = prob.core
	model = prob.model
	game_obj = prob.game_obj

	# Initialization
	prob.core.res .= 0.0
    stamp = VStamp()
	∇dyn = zeros(MMatrix{n,(n+m),T,n*(n+m)})

	# Cost
	cost_gradient!(game_obj, pdtraj)
	for i = 1:p
		# State cost
		for k = 1:N
			stampify!(stamp, :opt, i, :x, 1, k)
			n_obj = length(game_obj.E[i])
			for j = 1:n_obj
				valid(stamp, N, p) ? add2sub(core.res_sub[stamp], game_obj.E[i][j].cost[k].q) : nothing
			end
		end
		# Control Cost
		for k = 1:N-1
			stampify!(stamp, :opt, i, :u, i, k)
			n_obj = length(game_obj.E[i])
			for j = 1:n_obj
				valid(stamp, N, p) ? add2sub(core.res_sub[stamp], game_obj.E[i][j].cost[k].r[pu[i]]) : nothing
			end
		end
	end
	# Dynamics penalty
	for k = 1:N-1
		∇dynamics!(∇dyn, model, pdtraj, k)
		for i = 1:p
			λik = pdtraj.du[i][k]
			stampify!(stamp, :opt, i, :x, 1, k)
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], ∇dyn[:,core.dyn[:x][1]]'*λik) : nothing
			stampify!(stamp, :opt, i, :u, i, k)
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], ∇dyn[:,core.dyn[:u][i]]'*λik) : nothing
			stampify!(stamp, :opt, i, :x, 1, k+1)
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], -λik) : nothing
		end
	end

	# Constraints
	constraint_residual!(prob, pdtraj)

	# Dynamics
    for k = 1:N-1
        stampify!(stamp, :dyn, 1, :x, 1, k)
        valid(stamp, N, p) ? add2sub(core.res_sub[stamp], dynamics_residual(model, pdtraj, k)) : nothing # shouldn't allocate
    end
    return nothing
end

function regularize_residual!(core::NewtonCore, opts::Options, pdtraj::PrimalDualTraj, pdtraj_ref::PrimalDualTraj)
	N = core.probsize.N
	p = core.probsize.p
	pu = core.probsize.pu

	stamp = VStamp()
	for k = 1:N-1
		x = state(pdtraj.pr[k+1])
		x_ref = state(pdtraj_ref.pr[k+1])
		u = control(pdtraj.pr[k])
		u_ref = control(pdtraj_ref.pr[k])
		for i = 1:p
			stampify!(stamp, :opt, i, :x, 1, k+1)
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], opts.reg.x*(x - x_ref)) : nothing
			stampify!(stamp, :opt, i, :u, i, k)
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], opts.reg.u*(u[pu[i]] - u_ref[pu[i]])) : nothing
		end
	end
	return nothing
end


################################################################################
# Residual Jacobian
################################################################################

function residual_jacobian!(prob::GameProblem{KN,n,m,T,SVd,SVx}) where {KN,n,m,T,SVd,SVx}
	residual_jacobian!(prob, prob.pdtraj)
	return nothing
end

function residual_jacobian!(prob::GameProblem{KN,n,m,T,SVd,SVx},
	pdtraj::PrimalDualTraj{KN,n,m,T,SVd}) where {KN,n,m,T,SVd,SVx}
    N = prob.probsize.N
	p = prob.probsize.p
    pu = prob.probsize.pu
	model = prob.model
	core = prob.core
    game_obj = prob.game_obj

	# Reset!
	sparse_zero!(prob.core.jac)

    # Allocations
    stamp = Stamp()
    ∇dyn = zeros(MMatrix{n,(n+m),T,n*(n+m)})

	# Cost function
	cost_hessian!(game_obj, pdtraj)
	# Cost
	for i = 1:p
		# State cost
		for k = 1:N
			stampify!(stamp, :opt, i, :x, 1, k, :x, 1, k)
			n_obj = length(game_obj.E[i])
			for j = 1:n_obj
				valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], game_obj.E[i][j].cost[k].Q) : nothing
			end
		end
		# Control cost
		for k = 1:N-1
			stampify!(stamp, :opt, i, :u, i, k, :u, i, k)
			n_obj = length(game_obj.E[i])
			for j = 1:n_obj
				valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], game_obj.E[i][j].cost[k].R[pu[i],pu[i]]) : nothing
			end
		end
	end

	# Constraints
	constraint_jacobian_residual!(prob, pdtraj)

	# Dynamics
	for k = 1:N-1
        ∇dynamics!(∇dyn, model, pdtraj, k)
        # Bottom Left
        stampify!(stamp, :dyn, 1, :x, 1, k, :x, 1, k)
		∇dyn_x = ∇dyn[:,core.dyn[:x][1]]
		valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:x][1]]) : nothing
        for i = 1:p
            stampify!(stamp, :dyn, 1, :x, 1, k, :u, i, k)
            valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:u][i]]) : nothing
        end
        stampify!(stamp, :dyn, 1, :x, 1, k, :x, 1, k+1)
        valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], -1.0) : nothing
        # Top Right
        for i = 1:p
            stampify!(stamp, :opt, i, :x, 1, k, :λ, i, k)
            valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:x][1]]') : nothing
            stampify!(stamp, :opt, i, :u, i, k, :λ, i, k)
            valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:u][i]]') : nothing
            stampify!(stamp, :opt, i, :x, 1, k+1, :λ, i, k)
            valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], -1.0) : nothing
        end
    end
    return nothing
end

function regularize_residual_jacobian!(prob::GameProblem{KN,n,m,T,SVd,SVx}) where {KN,n,m,T,SVd,SVx}
	N = prob.probsize.N
	p = prob.probsize.p
	pu = prob.probsize.pu
	core = prob.core
	opts = prob.opts

	stamp = Stamp()
	for k = 1:N-1
		for i = 1:p
			stampify!(stamp, :opt, i, :x, 1, k+1, :x, 1, k+1)
			valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], opts.reg.x) : nothing
			stampify!(stamp, :opt, i, :u, i, k, :u, i, k)
			valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], opts.reg.u) : nothing
		end
	end
	return nothing
end


################################################################################
# Iterative Best Response Residual
################################################################################
function ibr_residual!(prob::GameProblem{KN,n,m,T,SVd,SVx}, i::Int) where {KN,n,m,T,SVd,SVx}
	ibr_residual!(prob, prob.pdtraj, i)
	return nothing
end

function ibr_residual!(prob::GameProblem{KN,n,m,T,SVd,SVx}, pdtraj::PrimalDualTraj{KN,n,m,T,SVd}, i::Int) where {KN,n,m,T,SVd,SVx}
	N = prob.probsize.N
	p = prob.probsize.p
    pu = prob.probsize.pu
	core = prob.core
	model = prob.model
	game_obj = prob.game_obj

	# Initialization
	prob.core.res .= 0.0
    stamp = VStamp()
	∇dyn = zeros(MMatrix{n,(n+m),T,n*(n+m)})

	# Cost
	cost_gradient!(game_obj, pdtraj, i)
	# State cost
	for k = 1:N
		stampify!(stamp, :opt, i, :x, 1, k)
		n_obj = length(game_obj.E[i])
		for j = 1:n_obj
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], game_obj.E[i][j].cost[k].q) : nothing
		end
	end
	# Control Cost
	for k = 1:N-1
		stampify!(stamp, :opt, i, :u, i, k)
		n_obj = length(game_obj.E[i])
		for j = 1:n_obj
			valid(stamp, N, p) ? add2sub(core.res_sub[stamp], game_obj.E[i][j].cost[k].r[pu[i]]) : nothing
		end
	end
	# Dynamics penalty
	for k = 1:N-1
		∇dynamics!(∇dyn, model, pdtraj, k)
		λik = pdtraj.du[i][k]
		stampify!(stamp, :opt, i, :x, 1, k)
		valid(stamp, N, p) ? add2sub(core.res_sub[stamp], ∇dyn[:,core.dyn[:x][1]]'*λik) : nothing
		stampify!(stamp, :opt, i, :u, i, k)
		valid(stamp, N, p) ? add2sub(core.res_sub[stamp], ∇dyn[:,core.dyn[:u][i]]'*λik) : nothing
		stampify!(stamp, :opt, i, :x, 1, k+1)
		valid(stamp, N, p) ? add2sub(core.res_sub[stamp], -λik) : nothing
	end

	# Constraints
	constraint_residual!(prob, pdtraj)

	# Dynamics
    for k = 1:N-1
        stampify!(stamp, :dyn, 1, :x, 1, k)
        valid(stamp, N, p) ? add2sub(core.res_sub[stamp], dynamics_residual(model, pdtraj, k)) : nothing # shouldn't allocate
    end
    return nothing
end

function ibr_regularize_residual!(core::NewtonCore, opts::Options, pdtraj::PrimalDualTraj, pdtraj_ref::PrimalDualTraj, i::Int)
	N = core.probsize.N
	p = core.probsize.p
	pu = core.probsize.pu

	stamp = VStamp()
	for k = 1:N-1
		x = state(pdtraj.pr[k+1])
		x_ref = state(pdtraj_ref.pr[k+1])
		u = control(pdtraj.pr[k])
		u_ref = control(pdtraj_ref.pr[k])
		stampify!(stamp, :opt, i, :x, 1, k+1)
		valid(stamp, N, p) ? add2sub(core.res_sub[stamp], opts.reg.x*(x - x_ref)) : nothing
		stampify!(stamp, :opt, i, :u, i, k)
		valid(stamp, N, p) ? add2sub(core.res_sub[stamp], opts.reg.u*(u[pu[i]] - u_ref[pu[i]])) : nothing
	end
	return nothing
end


################################################################################
# Iterative Best Response Residual Jacobian
################################################################################

function ibr_residual_jacobian!(prob::GameProblem{KN,n,m,T,SVd,SVx}, i::Int) where {KN,n,m,T,SVd,SVx}
	residual_jacobian!(prob, prob.pdtraj, i)
	return nothing
end

function ibr_residual_jacobian!(prob::GameProblem{KN,n,m,T,SVd,SVx},
	pdtraj::PrimalDualTraj{KN,n,m,T,SVd}, i::Int) where {KN,n,m,T,SVd,SVx}
    N = prob.probsize.N
	p = prob.probsize.p
    pu = prob.probsize.pu
	model = prob.model
	core = prob.core
    game_obj = prob.game_obj

	# Reset!
	sparse_zero!(prob.core.jac)

    # Allocations
    stamp = Stamp()
    ∇dyn = zeros(MMatrix{n,(n+m),T,n*(n+m)})

	# Cost function
	cost_hessian!(game_obj, pdtraj, i)
	# Cost
	# State cost
	for k = 1:N
		stampify!(stamp, :opt, i, :x, 1, k, :x, 1, k)
		n_obj = length(game_obj.E[i])
		for j = 1:n_obj
			valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], game_obj.E[i][j].cost[k].Q) : nothing
		end
	end
	# Control cost
	for k = 1:N-1
		stampify!(stamp, :opt, i, :u, i, k, :u, i, k)
		n_obj = length(game_obj.E[i])
		for j = 1:n_obj
			valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], game_obj.E[i][j].cost[k].R[pu[i],pu[i]]) : nothing
		end
	end

	# Constraints
	constraint_jacobian_residual!(prob, pdtraj)

	# Dynamics
	for k = 1:N-1
        ∇dynamics!(∇dyn, model, pdtraj, k)
        # Bottom Left
        stampify!(stamp, :dyn, 1, :x, 1, k, :x, 1, k)
		∇dyn_x = ∇dyn[:,core.dyn[:x][1]]
		valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:x][1]]) : nothing

        stampify!(stamp, :dyn, 1, :x, 1, k, :u, i, k)
        valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:u][i]]) : nothing

        stampify!(stamp, :dyn, 1, :x, 1, k, :x, 1, k+1)
        valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], -1.0) : nothing
        # Top Right
        stampify!(stamp, :opt, i, :x, 1, k, :λ, i, k)
        valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:x][1]]') : nothing
        stampify!(stamp, :opt, i, :u, i, k, :λ, i, k)
        valid(stamp, N, p) ? add2sub(core.jac_sub[stamp], ∇dyn[:,core.dyn[:u][i]]') : nothing
        stampify!(stamp, :opt, i, :x, 1, k+1, :λ, i, k)
        valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], -1.0) : nothing
    end
    return nothing
end

function ibr_regularize_residual_jacobian!(prob::GameProblem{KN,n,m,T,SVd,SVx}, i::Int) where {KN,n,m,T,SVd,SVx}
	N = prob.probsize.N
	p = prob.probsize.p
	pu = prob.probsize.pu
	core = prob.core
	opts = prob.opts

	stamp = Stamp()
	for k = 1:N-1
		stampify!(stamp, :opt, i, :x, 1, k+1, :x, 1, k+1)
		valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], opts.reg.x) : nothing
		stampify!(stamp, :opt, i, :u, i, k, :u, i, k)
		valid(stamp, N, p) ? addI2sub(core.jac_sub[stamp], opts.reg.u) : nothing
	end
	return nothing
end
