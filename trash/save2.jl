export
	rel_zinds,
	zinds

# For ∇c
function rel_zinds(con::ConstraintVals{T,Stage}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == 1
		return StaticArrays.SUnitRange(n+1,n+m) # n .+ 1:m
	elseif k == N
		return SOneTo(n) # 1:n
	else
		return SOneTo(n+m) # #[1:n; n .+ (1:m)]
	end
end

function rel_zinds(con::ConstraintVals{T,State}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == 1
		return @SVector zeros(0) # 1:0
	else
		return Zk._x # 1:n
	end
end

function rel_zinds(con::ConstraintVals{T,Control}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == N
		return @SVector zeros(0)
	else
		return SOneTo(m) # 1:m
	end
end

function rel_zinds(con::ConstraintVals{T,Coupled}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	n,m,N = size(solver)
	if k == 1
		return StaticArrays.SUnitRange(n+1,2*(n+m))
	elseif k == N-1
		return SOneTo(2*n+m) # [1:n; n .+ (1:m); n+m .+ (1:n)]
	elseif k == N
		return Zk._x # 1:n
	else
		return SOneTo(2*(n+m)) # [1:n; n .+ (1:m); n+m .+ (1:n); 2*n+m .+ (1:m)]
	end
end

function rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
	Zk::KnotPoint{T,N1,M,NM}, k::Int, n::Int, m::Int, pl::Vector{Vector{Int}}, p::Int, N::Int) where {T,Q<:TO.Implicit, N1, M, NM}
	if k == 1
		return StaticArrays.SUnitRnage(n+1,2*n+m) # [n .+ (1:m); n+m .+ (1:n)]
	elseif k == N
		return Zk._x # 1:n
	else
		return SOneTo(2*n+m) #[1:n; n .+ (1:m); n+m .+ (1:n)]
	end
end

# For ∇ci
function rel_zinds(con::ConstraintVals{T,Stage}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, i::Int, n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == 1
		return SVector{length(pu[i])}([n+j for j in pu[i]]) # n .+ pu[i]
	elseif k == N
		return Zk._x # 1:n
	else
		return SVector{n+length(pu[i])}([j for j in [1:n; n .+ pu[i]]]) # [1:n; n .+ pu[i]]
	end
end

function rel_zinds(con::ConstraintVals{T,State}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, i::Int, n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == 1
		return @SVector zeros(0) # 1:0
	else
		return Zk._x # 1:n
	end
end

function rel_zinds(con::ConstraintVals{T,Control}, Zk::KnotPoint{T,N1,M,NM},
	k::Int, i::Int, n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == N
		return @SVector zeros(0) # 1:0
	else
		return SVector{length(pu[1])}([j for j in pu[i]]) # pu[i]
	end
end

function rel_zinds(con::ConstraintVals{T,Coupled}, Zk::KnotPoint{T,N1,M,NM}, k::Int, i::Int, n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T,N1,M,NM}
	if k == 1
		return SVector{n+2*length(pu[i])}([j for j in [n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]]) # [n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
	elseif k == N-1
		return SVector{2*n+length(pu[i])}([j for j in [1:n; n .+ pu[i]; n+m .+ (1:n)]]) # [1:n; n .+ pu[i]; n+m .+ (1:n)]
	elseif k == N
		return SVector{n}([i for i=1:n]) # 1:n
	else
		return SVector{2*(n+length(pu[i]))}([j for j in [1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]]) # [1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]
	end
end

function rel_zinds(con::Union{ConstraintVals{T,Dynamical}, ConstraintVals{T,Coupled,<:DynamicsConstraint{Q}}},
	Zk::KnotPoint{T,N1,M,NM}, k::Int, i::Int, n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T,Q<:TO.Implicit, N1, M, NM}
	if k == 1
		return SVector{n+length(pu[i])}([j for j in [n .+ pu[i]; n+m .+ (1:n)]]) # [n .+ pu[i]; n+m .+ (1:n)]
	elseif k == N
		return SVector{n}([i for i=1:n]) # 1:n
	else
		return SVector{2*n+length(pu[i])}([j for j in [1:n; n .+ pu[i]; n+m .+ (1:n)]]) # [1:n; n .+ pu[i]; n+m .+ (1:n)]
	end
end

#
# static_base =
# StaticArrays.SUnitRange(n+1,n+m)
# SOneTo(n)
# SOneTo(n+m)
#
# @SVector zeros(0)
# SOneTo(n)
#
# @SVector zeros(0)
# SoneTo(m)
#
# StaticArrays.SUnitRange(n+1,2*(n+m))
# SOneTo(2*n+m)
# SOneTo(n)
# SOneTo(2*(n+m))
#
# StaticArrays(n+1,2*n+m)
# SOneTo(n)
# SOneTo(2*n+m)
#
#
#
# static_base_i =
# SVector{length(pu[i])}(n.+ pu[i])
# SOneTo(n)
# SVector{n+length(pu[i])}([1:n; n .+ pu[i]])
#
# @SVector zeros(0)
# SOneTo(n)
#
# @SVector zeros(0)
# SVector{length(pu[1])}(pu[i])
#
# SVector{n+2*length(pu[i])}([n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]])
# SVector{2*n+length(pu[i])}([1:n; n .+ pu[i]; n+m .+ (1:n)])
# SVector{n}([i for i=1:n])
# SVector{2*(n+length(pu[i]))}([1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]])
#
# SVector{n+length(pu[i])}([n .+ pu[i]; n+m .+ (1:n)])
# SVector{n}([i for i=1:n])
# SVector{2*n+length(pu[i])}([1:n; n .+ pu[i]; n+m .+ (1:n)])



struct StaticInds
	stage::Vector{AbstractVector{Int}}
	state::Vector{AbstractVector{Int}}
	control::Vector{AbstractVector{Int}}
	coupled::Vector{AbstractVector{Int}}
	dynamical::Vector{AbstractVector{Int}}
end


struct StaticIndsP
	stage::Vector{Vector{AbstractVector{Int}}}
	state::Vector{Vector{AbstractVector{Int}}}
	control::Vector{Vector{AbstractVector{Int}}}
	coupled::Vector{Vector{AbstractVector{Int}}}
	dynamical::Vector{Vector{AbstractVector{Int}}}
end

# sstage = [StaticArrays.SUnitRange(n+1,n+m), SOneTo(n), SOneTo(n+m)]
# sstate = [SVector{0}(zeros(0)), SOneTo(n)]
# scontrol = [SVector{0}(zeros(0)), SOneTo(m)]
# scoupled = [StaticArrays.SUnitRange(n+1,2*(n+m)), SOneTo(2*n+m), SOneTo(n), SOneTo(2*(n+m))]
# sdynamical = [StaticArrays.SUnitRange(n+1,2*n+m), SOneTo(n), SOneTo(2*n+m)]
#
# sstage_p = [SVector{length(pu[i])}(n.+ pu[i]), SOneTo(n), SVector{n+length(pu[i])}([1:n; n .+ pu[i]])]
# sstate_p = [SVector{0}(zeros(0)), SOneTo(n)]
# scontrol_p = [SVector{0}(zeros(0)), SVector{length(pu[1])}(pu[i])]
# scoupled_p = [
# 	SVector{n+2*length(pu[i])}([n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]]),
# 	SVector{2*n+length(pu[i])}([1:n; n .+ pu[i]; n+m .+ (1:n)]),
# 	SVector{n}([i for i=1:n]),
# 	SVector{2*(n+length(pu[i]))}([1:n; n .+ pu[i]; n+m .+ (1:n); 2*n+m .+ pu[i]])]
# sdynamical_p = [
# 	SVector{n+length(pu[i])}([n .+ pu[i]; n+m .+ (1:n)]),
# 	SVector{n}([i for i=1:n]),
# 	SVector{2*n+length(pu[i])}([1:n; n .+ pu[i]; n+m .+ (1:n)])]
#
#
# sinds = StaticInds5(sstage, sstate, scontrol, scoupled, sdynamical)
# sinds_p = StaticInds5(sstage_p, sstate_p, scontrol_p, scoupled_p, sdynamical_p)


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
