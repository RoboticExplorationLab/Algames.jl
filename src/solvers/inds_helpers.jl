export
	rel_zinds,
	zinds

# For ∇c
function rel_zinds(con::ConstraintVals{T,Stage}, sinds::StaticInds,
	k::Int, n::Int, m::Int, N::Int) where {T}
	if k == 1
		return n .+ sinds.sm # n .+ 1:m
	elseif k == N
		return sinds.sn # 1:n
	else
		return sinds.snm # #[1:n; n .+ (1:m)]
	end
end

function rel_zinds(con::ConstraintVals{T,State}, sinds::StaticInds,
	k::Int, n::Int, m::Int, N::Int) where {T}
	if k == 1
		return sinds.s0 # 1:0
	else
		return sinds.sn # 1:n
	end
end

function rel_zinds(con::ConstraintVals{T,Control}, sinds::StaticInds,
	k::Int, n::Int, m::Int, N::Int) where {T}
	if k == N
		return sinds.s0
	else
		return sinds.sm # 1:m
	end
end

function rel_zinds(con::ConstraintVals{T,Coupled},sinds::StaticInds,
	k::Int, n::Int, m::Int, N::Int) where {T}
	if k == 1
		return [n .+ sinds.sm; (n+m) .+ sinds.snm]
	elseif k == N-1
		return sinds.s2nm # [1:n; n .+ (1:m); n+m .+ (1:n)]
	elseif k == N
		return sinds.sn # 1:n
	else
		return sinds.s2n2m # [1:n; n .+ (1:m); n+m .+ (1:n); 2*n+m .+ (1:m)]
	end
end

function rel_zinds(con::Union{ConstraintVals{T,Dynamical},
	ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
	sinds::StaticInds, n::Int, m::Int, N::Int) where {T,Q<:TO.Implicit}
	if k == 1
		return n .+ sinds.snm # [n .+ (1:m); n+m .+ (1:n)]
	elseif k == N
		return sinds.sn # 1:n
	else
		return sinds.s2nm #[1:n; n .+ (1:m); n+m .+ (1:n)]
	end
end

# For ∇ci
function rel_zinds(con::ConstraintVals{T,Stage}, solver::DirectGamesSolver{T},
	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
	if k == 1
		return [n .+ solver.spu[i]] # n .+ pu[i]
	elseif k == N
		return solver.sinds.sn # 1:n
	else
		return [solver.sinds.sn; n.+ solver.spu[i]] # [1:n; n .+ pu[i]]
	end
end

function rel_zinds(con::ConstraintVals{T,State}, solver::DirectGamesSolver{T},
	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
	if k == 1
		return solver.sinds.s0 # 1:0
	else
		return solver.sinds.sn # 1:n
	end
end

function rel_zinds(con::ConstraintVals{T,Control}, solver::DirectGamesSolver{T},
	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
	if k == N
		return solver.sinds.s0 # 1:0
	else
		return solver.spu[i] # pu[i]
	end
end

function rel_zinds(con::ConstraintVals{T,Coupled}, solver::DirectGamesSolver{T},
	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
	if k == 1
		return  [n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn; (2*n+m) .+ solver.spu[i]] # [n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
	elseif k == N-1
		return  [solver.sinds.sn; n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn] # [1:n; n .+ pu[i]; n+m .+ (1:n)]
	elseif k == N
		return solver.sinds.sn # 1:n
	else
		return [solver.sinds.sn; n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn; (2*n+m) .+ solver.spu[i]] # [1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
	end
end

function rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
	solver::DirectGamesSolver{T}, k::Int, i::Int, n::Int, m::Int, N::Int) where {T,Q<:TO.Implicit}
	if k == 1
		return [n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn] # [n .+ pu[i]; n+m .+ (1:n)]
	elseif k == N
		return solver.sinds.sn # 1:n
	else
		return [solver.sinds.sn; n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn] # [1:n; n .+ pu[i]; n+m .+ (1:n)]
	end
end


# For g_
function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Stage}, k, i) where T
	n,m,N = size(solver)
	if k == 1
		return solver.uinds_p[i][k]
	elseif k == N
		return solver.xinds_p[i][k]
	else
		return [solver.xinds_p[i][k]; solver.uinds_p[i][k]]
	end
end

function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State}, k, i) where T
	n,m,N = size(solver)
	if k == 1
		return 1:0
	else
		return solver.xinds_p[i][k]
	end
end

function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Control}, k, i) where T
	n,m,N = size(solver)
	if k == N
		return 1:0
	else
		return solver.uinds_p[i][k]
	end
end

function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Coupled}, k, i) where T
	n,m,N = size(solver)
	if k == 1
		return [solver.uinds_p[i][k]; solver.xinds_p[i][k+1]; solver.uinds_p[i][k+1]]
	elseif k == N-1
		return [solver.xinds_p[i][k]; solver.uinds_p[i][k]; solver.xinds_p[i][k+1]]
	elseif k == N
		return [solver.xinds_p[i][k]]
	else
		return [solver.xinds_p[i][k]; solver.uinds_p[i][k]; solver.xinds_p[i][k+1]; solver.uinds_p[i][k+1]]
	end
end

function zinds(solver::DirectGamesSolver{T},
	con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
	k, i) where {T,Q<:TO.Implicit}
	n,m,N = size(solver)
	if k == 1
		return [solver.uinds_p[i][k]; solver.xinds_p[i][k+1]]
	elseif k == N
		return [solver.xinds_p[i][k]]
	else
		return [solver.xinds_p[i][k]; solver.uinds_p[i][k]; solver.xinds_p[i][k+1]]
	end
end

# For H_
function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Stage}, k) where T
	n,m,N = size(solver)
	if k == 1
		return solver.uinds[k]
	elseif k == N
		return solver.xinds[k]
	else
		return [solver.xinds[k]; solver.uinds[k]]
	end
end

function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State}, k) where T
	n,m,N = size(solver)
	if k == 1
		return 1:0
	else
		return solver.xinds[k]
	end
end

function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Control}, k) where T
	n,m,N = size(solver)
	if k == N
		return 1:0
	else
		return solver.uinds[k]
	end
end

function zinds(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Coupled}, k) where T
	n,m,N = size(solver)
	if k == 1
		return [solver.uinds[k]; solver.xinds[k+1]; solver.uinds[k+1]]
	elseif k == N-1
		return [solver.xinds[k]; solver.uinds[k]; solver.xinds[k+1]]
	elseif k == N
		return [solver.xinds[k]]
	else
		return [solver.xinds[k]; solver.uinds[k]; solver.xinds[k+1]; solver.uinds[k+1]]
	end
end

function zinds(solver::DirectGamesSolver{T},
	con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
	k) where {T,Q<:TO.Implicit}
	n,m,N = size(solver)
	if k == 1
		return [solver.uinds[k]; solver.xinds[k+1]]
	elseif k == N
		return [solver.xinds[k]]
	else
		return [solver.xinds[k]; solver.uinds[k]; solver.xinds[k+1]]
	end
end








#
# function rel_zinds(con::ConstraintVals{T,Stage}, sinds::StaticInds,
# 	k::Int, N::Int) where {T}
# 	if k == 1
# 		return sinds.stage[1] # n .+ 1:m
# 	elseif k == N
# 		return sinds.stage[2] # 1:n
# 	else
# 		return sinds.stage[3] # #[1:n; n .+ (1:m)]
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,State}, sinds::StaticInds,
# 	k::Int, N::Int) where {T}
# 	if k == 1
# 		return sinds.state[1] # 1:0
# 	else
# 		return sinds.state[2] # 1:n
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,Control}, sinds::StaticInds,
# 	k::Int, N::Int) where {T}
# 	if k == N
# 		return sinds.control[1]
# 	else
# 		return sinds.control[2] # 1:m
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,Coupled},sinds::StaticInds,
# 	k::Int, N::Int) where {T}
# 	if k == 1
# 		return sinds.coupled[1]
# 	elseif k == N-1
# 		return sinds.coupled[2] # [1:n; n .+ (1:m); n+m .+ (1:n)]
# 	elseif k == N
# 		return sinds.coupled[3] # 1:n
# 	else
# 		return sinds.coupled[4] # [1:n; n .+ (1:m); n+m .+ (1:n); 2*n+m .+ (1:m)]
# 	end
# end
#
# function rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
# 	sinds::StaticInds, N::Int) where {T,Q<:TO.Implicit}
# 	if k == 1
# 		return sinds.dynamical[1] # [n .+ (1:m); n+m .+ (1:n)]
# 	elseif k == N
# 		return sinds.dynamical[2] # 1:n
# 	else
# 		return sinds.dynamical[3] #[1:n; n .+ (1:m); n+m .+ (1:n)]
# 	end
# end


# function rel_zinds(con::ConstraintVals{T,State}, Zk::KnotPoint{T,N1,M,NM},
# 	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
# 	if k == 1
# 		return @SVector zeros(0) # 1:0
# 	else
# 		return Zk._x # 1:n
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,Control}, Zk::KnotPoint{T,N1,M,NM},
# 	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
# 	if k == N
# 		return @SVector zeros(0)
# 	else
# 		return SOneTo(m) # 1:m
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,Coupled}, Zk::KnotPoint{T,N1,M,NM},
# 	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
# 	n,m,N = size(solver)
# 	if k == 1
# 		return StaticArrays.SUnitRange(n+1,2*(n+m))
# 	elseif k == N-1
# 		return SOneTo(2*n+m) # [1:n; n .+ (1:m); n+m .+ (1:n)]
# 	elseif k == N
# 		return Zk._x # 1:n
# 	else
# 		return SOneTo(2*(n+m)) # [1:n; n .+ (1:m); n+m .+ (1:n); 2*n+m .+ (1:m)]
# 	end
# end
#
# function rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
# 	Zk::KnotPoint{T,N1,M,NM}, k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,Q<:TO.Implicit, N1, M, NM}
# 	if k == 1
# 		return StaticArrays.SUnitRnage(n+1,2*n+m) # [n .+ (1:m); n+m .+ (1:n)]
# 	elseif k == N
# 		return Zk._x # 1:n
# 	else
# 		return SOneTo(2*n+m) #[1:n; n .+ (1:m); n+m .+ (1:n)]
# 	end
# end
#
#
# inds_stage = [
# 	StaticArrays.SUnitRange(n+1,n+m),
# 	SOneTo(n),
# 	SOneTo(n+m)]
# inds_state = [
# 	SVector{0}(zeros(0)),
# 	SOneTo(n)]
# inds_control = [
# 	SVector{0}(zeros(0)),
# 	SOneTo(m)]
# inds_coupled = [
# 	StaticArrays.SUnitRange(n+1,2*(n+m)),
# 	SOneTo(2*n+m),
# 	SOneTo(n),
# 	SOneTo(2*(n+m))]
# inds_dynamical = [
# 	StaticArrays.SUnitRange(n+1,2*n+m),
# 	SOneTo(n),
# 	SOneTo(2*n+m)]
#
# inds_stage_p = [[
# 	SVector{length(pu[i])}(n.+ pu[i]),
# 	SOneTo(n),
# 	SVector{n+length(pu[i])}([1:n; n .+ pu[i]])] for i=1:p]
# inds_state_p = [[
# 	SVector{0}(zeros(0)),
# 	SOneTo(n)] for i=1:p]
# inds_control_p = [[
# 	SVector{0}(zeros(0)),
# 	SVector{length(pu[1])}(pu[i])] for i=1:p]
# inds_coupled_p = [[
# 	SVector{n+2*length(pu[i])}([n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]),
# 	SVector{2*n+length(pu[i])}([1:n; n .+ pu[i]; n+m .+ (1:n)]),
# 	SVector{n}([i for i=1:n]),
# 	SVector{2*(n+length(pu[i]))}([1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]])] for i=1:p]
# inds_dynamical_p = [[
# 	SVector{n+length(pu[i])}([n .+ pu[i]; n+m .+ (1:n)]),
# 	SVector{n}([i for i=1:n]),
# 	SVector{2*n+length(pu[i])}([1:n; n .+ pu[i]; n+m .+ (1:n)])] for i=1:p]
# sinds = StaticInds(inds_stage, inds_state, inds_control, inds_coupled, inds_dynamical)
# sinds_p = StaticIndsP(inds_stage_p, inds_state_p, inds_control_p, inds_coupled_p, inds_dynamical_p)





# For ∇ci
# function con_rel_zinds(con::ConstraintVals{T,Stage}, solver::DirectGamesSolver{T},
# 	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == 1
# 		@show size(con.∇c[k])
# 		return con.∇c[k][:, n .+ solver.spu[i]] # n .+ pu[i]
# 	elseif k == N
# 		return con.∇c[k][:, solver.sinds.sn] # 1:n
# 	else
# 		return con.∇c[k][:, [solver.sinds.sn; n.+ solver.spu[i]]] # [1:n; n .+ pu[i]]
# 	end
# end
#
# function con_rel_zinds(con::ConstraintVals{T,State}, solver::DirectGamesSolver{T},
# 	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == 1
# 		return con.∇c[k][:, solver.sinds.s0] # 1:0
# 	else
# 		return con.∇c[k][:, solver.sinds.sn] # 1:n
# 	end
# end
#
# function con_rel_zinds(con::ConstraintVals{T,Control}, solver::DirectGamesSolver{T},
# 	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == N
# 		return con.∇c[k][:, solver.sinds.s0] # 1:0
# 	else
# 		return con.∇c[k][:, solver.spu[i]] # pu[i]
# 	end
# end
#
# function con_rel_zinds(con::ConstraintVals{T,Coupled}, solver::DirectGamesSolver{T},
# 	k::Int, i::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == 1
# 		return  con.∇c[k][:, [n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn; (2*n+m) .+ solver.spu[i]]] # [n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
# 	elseif k == N-1
# 		return  con.∇c[k][:, [solver.sinds.sn; n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn]] # [1:n; n .+ pu[i]; n+m .+ (1:n)]
# 	elseif k == N
# 		return con.∇c[k][:, solver.sinds.sn] # 1:n
# 	else
# 		return con.∇c[k][:, [solver.sinds.sn; n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn; (2*n+m) .+ solver.spu[i]]] # [1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
# 	end
# end
#
# function con_rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
# 	solver::DirectGamesSolver{T}, k::Int, i::Int, n::Int, m::Int, N::Int) where {T,Q<:TO.Implicit}
# 	if k == 1
# 		return con.∇c[k][:, [n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn]] # [n .+ pu[i]; n+m .+ (1:n)]
# 	elseif k == N
# 		return con.∇c[k][:, solver.sinds.sn] # 1:n
# 	else
# 		return con.∇c[k][:, [solver.sinds.sn; n .+ solver.spu[i]; (n+m) .+ solver.sinds.sn]] # [1:n; n .+ pu[i]; n+m .+ (1:n)]
# 	end
# end



# struct StaticInds
# 	stage::Vector{AbstractVector{Int}}
# 	state::Vector{AbstractVector{Int}}
# 	control::Vector{AbstractVector{Int}}
# 	coupled::Vector{AbstractVector{Int}}
# 	dynamical::Vector{AbstractVector{Int}}
# end
#
# struct StaticIndsP
# 	stage::Vector{Vector{AbstractVector{Int}}}
# 	state::Vector{Vector{AbstractVector{Int}}}
# 	control::Vector{Vector{AbstractVector{Int}}}
# 	coupled::Vector{Vector{AbstractVector{Int}}}
# 	dynamical::Vector{Vector{AbstractVector{Int}}}
# end

# For ∇c
# function con_rel_zinds(con::ConstraintVals{T,Stage}, sinds::StaticInds,
# 	k::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == 1
# 		return con.∇c[k][:, n .+ sinds.sm] # n .+ 1:m
# 	elseif k == N
# 		return con.∇c[k][:, sinds.sn] # 1:n
# 	else
# 		return con.∇c[k][:, sinds.snm] # [1:n; n .+ (1:m)]
# 	end
# end
#
# function con_rel_zinds(con::ConstraintVals{T,State}, sinds::StaticInds,
# 	k::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == 1
# 		return con.∇c[k][:, sinds.s0] # 1:0
# 	else
# 		return con.∇c[k][:, sinds.sn] # 1:n
# 	end
# end
#
# function con_rel_zinds(con::ConstraintVals{T,Control}, sinds::StaticInds,
# 	k::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == N
# 		return con.∇c[k][:, sinds.s0]
# 	else
# 		return con.∇c[k][:, sinds.sm] # 1:m
# 	end
# end
#
# function con_rel_zinds(con::ConstraintVals{T,Coupled},sinds::StaticInds,
# 	k::Int, n::Int, m::Int, N::Int) where {T}
# 	if k == 1
# 		return con.∇c[k][:, [n .+ sinds.sm; (n+m) .+ sinds.snm]]
# 	elseif k == N-1
# 		return con.∇c[k][:, sinds.s2nm] # [1:n; n .+ (1:m); n+m .+ (1:n)]
# 	elseif k == N
# 		return con.∇c[k][:, sinds.sn] # 1:n
# 	else
# 		return con.∇c[k][:, sinds.s2n2m] # [1:n; n .+ (1:m); n+m .+ (1:n); 2*n+m .+ (1:m)]
# 	end
# end
#
# function con_rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
# 	sinds::StaticInds, n::Int, m::Int, N::Int) where {T,Q<:TO.Implicit}
# 	if k == 1
# 		return con.∇c[k][:, n .+ sinds.snm] # [n .+ (1:m); n+m .+ (1:n)]
# 	elseif k == N
# 		return con.∇c[k][:, sinds.sn] # 1:n
# 	else
# 		return con.∇c[k][:, sinds.s2nm] #[1:n; n .+ (1:m); n+m .+ (1:n)]
# 	end
# end


# k = N
# i = 2
# con = solver_directgames.constraints.constraints[end]
# con = solver_directgames.constraints.constraints[1]
# @allocated con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @btime con_rel_zinds($con,$solver_directgames.sinds,$k,$n,$m,$N)
# test_allocation(con_rel_zinds, (con,solver_directgames.sinds,k,n,m,N))
#
#
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = rel_zinds(con,solver_directgames.sinds,k,n,m,N)
# @allocated rel_zind_i = rel_zinds(con,solver_directgames.sinds,k,n,m,N)












#
# # For ∇ci
# function rel_zinds(con::ConstraintVals{T,Stage}, sinds_p::StaticIndsP,
# 	k::Int, i::Int, N::Int) where {T}
# 	if k == 1
# 		return sinds_p.stage[i][1] # n .+ pu[i]
# 	elseif k == N
# 		return sinds_p.stage[i][2] # 1:n
# 	else
# 		return sinds_p.stage[i][3] # [1:n; n .+ pu[i]]
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,State}, sinds_p::StaticIndsP,
# 	k::Int, i::Int, N::Int) where {T}
# 	if k == 1
# 		return sinds_p.state[i][1] # 1:0
# 	else
# 		return sinds_p.state[i][2] # 1:n
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,Control}, sinds_p::StaticIndsP,
# 	k::Int, i::Int, N::Int) where {T}
# 	if k == N
# 		return sinds_p.control[i][1] # 1:0
# 	else
# 		return sinds_p.control[i][2] # pu[i]
# 	end
# end
#
# function rel_zinds(con::ConstraintVals{T,Coupled}, sinds_p::StaticIndsP,
# 	k::Int, i::Int, N::Int) where {T}
# 	if k == 1
# 		return sinds_p.coupled[i][1] # [n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
# 	elseif k == N-1
# 		return sinds_p.coupled[i][2] # [1:n; n .+ pu[i]; n+m .+ (1:n)]
# 	elseif k == N
# 		return sinds_p.coupled[i][3] # 1:n
# 	else
# 		return sinds_p.coupled[i][4] # [1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
# 	end
# end
#
# function rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
# 	sinds_p::StaticIndsP, k::Int, i::Int, N::Int) where {T,Q<:TO.Implicit}
# 	if k == 1
# 		return sinds_p.dynamical[i][1] # [n .+ pu[i]; n+m .+ (1:n)]
# 	elseif k == N
# 		return sinds_p.dynamical[i][2] # 1:n
# 	else
# 		return sinds_p.dynamical[i][3] # [1:n; n .+ pu[i]; n+m .+ (1:n)]
# 	end
# end
