################################################################################
# PrimalDualTraj
################################################################################

mutable struct PrimalDualTraj{KN,n,m,T,SVd}
    probsize::ProblemSize
    pr::Traj{n,m,T,KN} # Primal trajectory
    du::Vector{SVd} # Dual trajectory
end

function PrimalDualTraj(probsize::ProblemSize, dt::T; f=rand, amplitude=1e-8) where {T}
    N = probsize.N
    n = probsize.n
    m = probsize.m
    p = probsize.p
    pr = Traj(n,m,dt,N)
    du = [[amplitude*f(SVector{n,T}) for k=1:N-1] for i=1:p]
    for k = 1:N
        pr[k].z = amplitude*f(n+m)
    end
    TYPE = (eltype(pr),n,m,T,eltype(du))
    return PrimalDualTraj{TYPE...}(probsize, pr, du)
end

################################################################################
# Methods
################################################################################

function init_traj!(pdtraj::PrimalDualTraj{KN,n,m,T,SVd}; x0=1e-8*rand(SVector{n,T}),
    s::Int=2^10, f=rand, amplitude=1e-8) where {KN,n,m,T,SVd}
    N = pdtraj.probsize.N
    p = pdtraj.probsize.p

    for k = 1:N
        pdtraj.pr[k].z = (k+s<=N) ? pdtraj.pr[k+s].z : amplitude*f(SVector{n+m,T})
    end
    for i = 1:p
        for k = 1:N-1
            pdtraj.du[i][k] = (k+s<=N-1) ? pdtraj.du[i][k+s] : amplitude*f(SVector{n,T})
        end
    end
    RobotDynamics.set_state!(pdtraj.pr[1], x0)
    return nothing
end

function set_traj!(core::NewtonCore, Δpdtraj::PrimalDualTraj{KN,n,m,T,SVd},
    Δtraj::AbstractVector) where {KN,n,m,T,SVd}
    N = core.probsize.N
    p = core.probsize.p
    pu = core.probsize.pu
    # Primals
    indu = zeros(Int,m)
	stamp = HStamp()
    for k = 1:N-1
        # States
		stampify!(stamp, :x, 1, k+1)
        ind = horizontal_idx(core, stamp)
        RobotDynamics.set_state!(Δpdtraj.pr[k+1], Δtraj[ind])
        # Controls
        for i = 1:p
			stampify!(stamp, :u, i, k)
            indu[pu[i]] = horizontal_idx(core, stamp)
        end
        RobotDynamics.set_control!(Δpdtraj.pr[k], Δtraj[indu])
    end
    # Duals
    for i = 1:p
        for k = 1:N-1
			stampify!(stamp, :λ, i, k)
            ind = horizontal_idx(core, stamp)
            Δpdtraj.du[i][k] = Δtraj[ind]
        end
    end
    return nothing
end


function get_traj!(core::NewtonCore, Δpdtraj::PrimalDualTraj{KN,n,m,T,SVd},
    Δtraj::AbstractVector) where {KN,n,m,T,SVd}
    N = core.probsize.N
    p = core.probsize.p
    pu = core.probsize.pu
    # Primals
	indu = zeros(Int,m)
	stamp = HStamp()
	for k = 1:N-1
		# States
		stampify!(stamp, :x, 1, k+1)
		ind = horizontal_idx(core, stamp)
		Δtraj[ind] = state(Δpdtraj.pr[k+1])
        # Controls
        for i = 1:p
			stampify!(stamp, :u, i, k)
            indu[pu[i]] = horizontal_idx(core, stamp)
        end
        Δtraj[indu] = control(Δpdtraj.pr[k])
    end
    # Duals
    for i = 1:p
        for k = 1:N-1
			stampify!(stamp, :λ, i, k)
            ind = horizontal_idx(core, stamp)
            Δtraj[ind] = Δpdtraj.du[i][k]
        end
    end
    return nothing
end

function update_traj!(target::PrimalDualTraj{KN,n,m,T,SVd},
    source::PrimalDualTraj{KN,n,m,T,SVd}, α::T,
    Δ::PrimalDualTraj{KN,n,m,T,SVd}) where {KN,n,m,T,SVd,SVx}
    N = target.probsize.N
    p = target.probsize.p
    # Primals
    for k = 1:N-1
        # States
        RobotDynamics.set_state!(target.pr[k+1], state(source.pr[k+1]) + α*state(Δ.pr[k+1]))
        # Controls
        RobotDynamics.set_control!(target.pr[k], control(source.pr[k]) + α*control(Δ.pr[k]))
    end
    # Duals
    for i = 1:p
        for k = 1:N-1
            target.du[i][k] = source.du[i][k] + α*Δ.du[i][k]
        end
    end
    return nothing
end

function Δ_step(pdtraj::PrimalDualTraj, α::T) where {T}
	s = 0.0
	N = pdtraj.probsize.N
	n = pdtraj.probsize.n
	m = pdtraj.probsize.m
	p = pdtraj.probsize.p
	for k = 1:N-1
		s += norm(state(pdtraj.pr[k+1]), 1) # xk+1
		s += norm(control(pdtraj.pr[k]), 1) # uk
		for i = 1:p
			# s += norm(pdtraj.du[i][k], 1) # λik
		end
	end
	s *= α
	# s /= (N-1)*(n+m+n*p)
	s /= (N-1)*(n+m)
	return s
end
