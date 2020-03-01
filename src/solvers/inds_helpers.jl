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
