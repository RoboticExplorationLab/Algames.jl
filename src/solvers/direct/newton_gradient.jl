export
	update_g_!,
	update_g_con!


# function update_g_old!(solver::DirectGamesSolver; use_copy=false::Bool)
# 	n,m,N = size(solver)
# 	n,m,pu,p = size(solver.model)
# 	if use_copy
# 		Z = solver.Z̄
# 		ν = solver.ν_
# 	else
# 		Z = solver.Z
# 		ν = solver.ν
# 	end
#
# 	g_ = solver.g_
# 	∇F = solver.∇F
# 	C = solver.C
# 	dyn = solver.dyn_constraints.constraints[1].vals
#
# 	# Fill g_
# 	for i = 1:p
# 		# Dynamics and cost
# 		# uk
# 		for k = 1:N-1
# 			ix = Z[k]._x
# 			iu = Z[k]._u
# 			fdu = ∇F[k][ix,iu][:,pu[i]] # Bki [n,mi]
# 			g_[solver.uinds_p[i][k]] = C[i].u[k] + fdu'*ν[i][k] #ok
# 		end
# 		# xk
# 		for k = 2:N-1
# 			ix = Z[k]._x
# 			fdx = ∇F[k][ix,ix] # Ak [n,n]
# 			g_[solver.xinds_p[i][k]] = C[i].x[k] + fdx'*ν[i][k] - ν[i][k-1] #ok
# 		end
# 		g_[solver.xinds_p[i][N]] = C[i].x[N] - ν[i][N-1] #ok
#
# 		# Constraints
# 	    for c in eachindex(solver.constraints.constraints)
# 			con = solver.constraints.constraints[c]
# 		# for con in solver.constraints.constraints
# 			if typeof(con).parameters[2] == State
# 				for k in intersect(con.inds,2:N)
# 					rel_zind_i = rel_zinds(con,solver.sinds_p,k,i,N)
# 					zind_i = zinds(solver,con,k,i)
# 					∇ci = con.∇c[k][:,rel_zind_i]
# 					Iμ = TO.penalty_matrix(con,k)
# 					@inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] # best
# 					# # view_g = solver.state_view_g[i][k-1] #bad
# 					# # view_g .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #bad
# 					# # solver.state_view_g[i][k-1] .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #Bad
# 				end
# 			else
# 				for k in con.inds
# 					rel_zind_i = rel_zinds(con,solver.sinds_p,k,i,N)
# 					zind_i = zinds(solver,con,k,i)
# 					∇ci = con.∇c[k][:,rel_zind_i]
# 					Iμ = TO.penalty_matrix(con,k)
# 					@inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k]
# 				end
# 			end
# 		end
#
# 		# Dynamics constraints
# 		for k = 1:N-1
# 			Iγi = Diagonal(solver.γ[i][k])
# 			ix = Z[k]._x
# 			iu = Z[k]._u
# 			fdui = ∇F[k][ix,iu][:,pu[i]] # Bki [n,mi]
# 			g_[solver.uinds_p[i][k]] += fdui'*Iγi*dyn[k]# need check
# 		end
# 		for k = 2:N-1 ##### MAYBE K = 3:N-1
# 			Iγi = Diagonal(solver.γ[i][k])
# 			ix = Z[k]._x
# 			fdx = ∇F[k][ix,ix] # Ak [n,n]
# 			g_[solver.xinds_p[i][k]] += fdx'*Iγi*dyn[k]# need check
# 		end
# 		for k = 1:N-1
# 			Iγi = Diagonal(solver.γ[i][k])
# 			g_[solver.xinds_p[i][k+1]] += -Iγi*dyn[k]# need check
# 		end
# 	end
# 	# Dynamics Constraints
# 	col_off = (n*p+m)*(N-1) #ok
# 	for k = 1:N-1
# 		g_[col_off+(k-1)*n .+ (1:n)] = dyn[k] #ok
# 	end
# 	return nothing
# end

function update_g_!(solver::DirectGamesSolver; use_copy=false::Bool)
	dyn = solver.dyn_constraints.constraints[1].vals
	solver.g_ .*= 0.0 #####################################################################################
	update_g_!(solver,dyn)
	return nothing
end

function update_g_!(solver::DirectGamesSolver, dyn::Vector{SVector{nn,T}}; use_copy=false::Bool) where {nn,T}
	n,m,N = size(solver)
	n,m,pu,p = size(solver.model)
	if use_copy
		Z = solver.Z̄
		ν = solver.ν_
	else
		Z = solver.Z
		ν = solver.ν
	end

	g_ = solver.g_
	∇F = solver.∇F
	C = solver.C

	# Fill g_
	for i = 1:p
		# Dynamics and cost
		# uk
		spui = solver.spu[i]
		for k = 1:N-1
			ix = Z[k]._x
			iu = Z[k]._u
			fdu = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
			g_[solver.uinds_p[i][k]] = C[i].u[k] + fdu'*ν[i][k] #ok
		end
		# xk
		for k = 2:N-1
			ix = Z[k]._x
			fdx = ∇F[k][ix,ix] # Ak [n,n]
			g_[solver.xinds_p[i][k]] = C[i].x[k] + fdx'*ν[i][k] - ν[i][k-1] #ok
		end
		g_[solver.xinds_p[i][N]] = C[i].x[N] - ν[i][N-1] #ok
		########################################################################
		# MAYBE USEFUL
		########################################################################
		# # Dynamics constraints
		# for k = 1:N-1
		# 	Iγi = Diagonal(solver.γ[i][k])
		# 	ix = Z[k]._x
		# 	iu = Z[k]._u
		# 	fdui = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
		# 	g_[solver.uinds_p[i][k]] += fdui'*Iγi*dyn[k]# need check
		# end
		# for k = 2:N-1 ##### MAYBE K = 3:N-1
		# 	Iγi = Diagonal(solver.γ[i][k])
		# 	ix = Z[k]._x
		# 	fdx = ∇F[k][ix,ix] # Ak [n,n]
		# 	g_[solver.xinds_p[i][k]] += fdx'*Iγi*dyn[k]# need check
		# end
		# for k = 1:N-1
		# 	Iγi = Diagonal(solver.γ[i][k])
		# 	g_[solver.xinds_p[i][k+1]] += -Iγi*dyn[k]# need check
		# end
		########################################################################
		########################################################################
		########################################################################
	end
	for c in eachindex(solver.constraints.constraints)
		con = solver.constraints.constraints[c]
		update_g_con!(solver,con,n,m,p,N)
	end

	# Dynamics Constraints
	col_off = (n*p+m)*(N-1) #ok
	for k = 1:N-1
		g_[col_off+(k-1)*n .+ (1:n)] = dyn[k] #ok
	end
	return nothing
end

function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Stage},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	xinds_p = solver.xinds_p
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k == 1
				# rel_zind_i = [n .+ spu[i]]
				rel_zind_i = n .+ spu[i]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = uinds_p[i][k]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			elseif k == N
				rel_zind_i = sinds.sn
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = xinds_p[i][k]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			else
				rel_zind_i = [sinds.sn; n.+ spu[i]]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [xinds_p[i][k]; uinds_p[i][k]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			end
		end
	end
	return nothing
end

function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State},
	n::Int, m::Int, p::Int, N::Int) where {T}
	xinds_p = solver.xinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k >= 2
				rel_zind_i = sinds.sn
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = xinds_p[i][k]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			end
		end
	end
	return nothing
end

function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Control},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k <= N-1
				rel_zind_i = spu[i]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = uinds_p[i][k]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			end
		end
	end
	return nothing
end

function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Coupled},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	xinds_p = solver.xinds_p
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k == 1
				rel_zind_i = [n .+ spu[i]; (n+m) .+ sinds.sn; (2*n+m) .+ spu[i]]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [uinds_p[i][k]; xinds_p[i][k+1]; uinds_p[i][k+1]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			elseif k == N-1
				rel_zind_i = [sinds.sn; n .+ spu[i]; (n+m) .+ sinds.sn]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			elseif k == N
				rel_zind_i = sinds.sn
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [xinds_p[i][k]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			else
				rel_zind_i = [sinds.sn; n .+ spu[i]; (n+m) .+ sinds.sn; (2*n+m) .+ spu[i]]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]; uinds_p[i][k+1]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			end
		end
	end
	return nothing
end

function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Dynamical},
	n::Int, m::Int, p::Int, N::Int) where {T}
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k == 1
				rel_zind_i = [n .+ spu[i]; (n+m) .+ sinds.sn]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [uinds_p[i][k]; xinds_p[i][k+1]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			elseif k == N
				rel_zind_i = sinds.sn
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [xinds_p[i][k]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			else
				rel_zind_i = [sinds.sn; n .+ spu[i]; (n+m) .+ sinds.sn]
				∇ci = con.∇c[j][:,rel_zind_i]
				zind_i = [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]]
				Iμ = TO.penalty_matrix(con,j)
				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
			end
		end
	end
	return nothing
end
