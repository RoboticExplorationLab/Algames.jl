export
	update_H_!,
	update_H_con!,
	set_view_H!


# function update_H_old!(solver::DirectGamesSolver)
# 	n,m,N = size(solver)
# 	n,m,pu,p = size(solver.model)
#
# 	Z = solver.Z
# 	H_ = solver.H_
# 	∇F = solver.∇F
# 	C = solver.C
#
# 	# Fill H_
# 	# Upper-left block
# 	for i = 1:p
# 		H_[solver.uinds_p[i][1], solver.uinds_H[i][1]] .= C[i].uu[1] #ok
# 		for k = 2:N-1
# 			H_[solver.xinds_p[i][k], solver.xinds[k]] .= C[i].xx[k] #ok
# 			H_[solver.uinds_p[i][k], solver.xinds[k]] .= C[i].ux[k] #ok
# 			H_[solver.xinds_p[i][k], solver.uinds_H[i][k]] .= C[i].ux[k]' #ok
# 			H_[solver.uinds_p[i][k], solver.uinds_H[i][k]] .= C[i].uu[k] #ok
# 		end
# 		H_[solver.xinds_p[i][N], solver.xinds[N]] .= C[i].xx[N] #ok
#
# 		for c in eachindex(solver.constraints.constraints)
# 			con = solver.constraints.constraints[c]
# 		# for con in solver.constraints.constraints
# 			if typeof(con).parameters[2] == State
# 				for (j,k) in enumerate(intersect(con.inds,2:N))
# 					rel_zind_i = rel_zinds(con,solver,k,i,n,m,N)
# 					rel_zind = rel_zinds(con,solver.sinds,k,n,m,N)
# 					∇ci = con.∇c[j][:,rel_zind_i] #ok
# 					∇c = con.∇c[j][:,rel_zind] #ok
# 					Iμ = TO.penalty_matrix(con,j)
# 					# H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
# 					# solver.state_view_H[i][k-1] .+= Sparse(∇ci'*Iμ*∇c) #
# 					solver.state_view_H[i][k-1] .+= ∇ci'*Iμ*∇c # Best
# 				end
# 			else
# 				for (j,k) in enumerate(con.inds)
# 					rel_zind_i = rel_zinds(con,solver,k,i,n,m,N)
# 					rel_zind = rel_zinds(con,solver.sinds,k,n,m,N)
# 					zind_i = zinds(solver,con,k,i)
# 					zind = zinds(solver,con,k)
# 					∇ci = con.∇c[j][:,rel_zind_i] #ok
# 					∇c = con.∇c[j][:,rel_zind] #ok
# 					Iμ = TO.penalty_matrix(con,j)
# 					H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
# 				end
# 			end
#         end
# 	end
#
#
# 		# # Dynamics Constraints
# 		# for k = 1:N-1
# 		# 	ix = Z[k]._x
# 		# 	iu = Z[k]._u
# 		# 	fdx = ∇F[k][ix,ix] # Ak [n,n]
# 		# 	fdu = ∇F[k][ix,iu] # Bk [n,m]
# 		# 	fdui = ∇F[k][ix,iu][:,pu[i]] # Bki [n,mi]
# 		# 	Iγi = Diagonal(solver.γ[i][k])
# 		# 	if k >= 2
# 		# 		H_[solver.xinds_p[i][k]  , solver.xinds[k]  ] .+=  fdx' *Iγi*fdx
# 		# 		H_[solver.xinds_p[i][k]  , solver.uinds[k]  ] .+=  fdx' *Iγi*fdu
# 		# 		H_[solver.xinds_p[i][k]  , solver.xinds[k+1]] .+= -fdx' *Iγi
# 		#
# 		# 		H_[solver.uinds_p[i][k]  , solver.xinds[k]  ] .+=  fdui'*Iγi*fdx
# 		#
# 		# 		H_[solver.xinds_p[i][k+1], solver.xinds[k]  ] .+= -      Iγi*fdx
# 		# 	end
# 		# 	H_[solver.uinds_p[i][k]  , solver.uinds[k]  ] .+=  fdui'*Iγi*fdu
# 		# 	H_[solver.uinds_p[i][k]  , solver.xinds[k+1]] .+= -fdui'*Iγi
# 		#
# 		# 	H_[solver.xinds_p[i][k+1], solver.uinds[k]  ] .+= -      Iγi*fdu
# 		# 	H_[solver.xinds_p[i][k+1], solver.xinds[k+1]] .+=        Iγi
# 		# end
#
#
# 		# for con in solver.dyn_constraints.constraints
# 		# 	for k in con.inds
# 		# 		rel_zind_i = rel_zinds(con,k,i,n,m,pu,p,N)
# 		# 		rel_zind = rel_zinds(con,k,n,m,pu,p,N)
# 		# 		zind_i = zinds(solver,con,k,i)
# 		# 		zind = zinds(solver,con,k)
# 		# 		∇c = con.∇c[k]
# 		# 		∇ci = con.∇c[k][:,rel_zind]
# 		# 		Iγi = Diagonal(solver.γ[i][k])
# 		# 		H_[zind_i,zind] += ∇ci'*Iγi*∇c# need check
#         #     end
#         # end ######## need to check all indices
# 	#
# 	# end
#
# 	# Upper-right block
# 	for i = 1:p
# 		# pli = @SVector [j for j in pl[i]]
# 		for k = 2:N-1
# 			ix = Z[k]._x #ok
# 			fdx = ∇F[k][ix,ix] #ok
# 			H_[solver.xinds_p[i][k], solver.νinds[i][k]] .= fdx' #ok
# 			H_[solver.xinds_p[i][k], solver.νinds[i][k-1]] .= -Diagonal(@SVector ones(n)) #ok
# 		end
# 		H_[solver.xinds_p[i][N], solver.νinds[i][N-1]] .= -Diagonal(@SVector ones(n)) #ok
#
# 		for k = 1:N-1
# 			ix = Z[k]._x #ok
# 			iu = Z[k]._u #ok
# 			fdu = ∇F[k][ix,iu][:,pu[i]] #ok
# 			H_[solver.uinds_p[i][k], solver.νinds[i][k]] .= fdu' #ok
# 		end
# 	end
#
# 	# Lower-left block
# 	for k = 2:N-1
# 		ix = Z[k]._x #ok
# 		fdx = ∇F[k][ix,ix] #ok
# 		H_[solver.νinds_p[k-1], solver.xinds[k]] .= -Diagonal(@SVector ones(n)) #ok
# 		H_[solver.νinds_p[k], solver.xinds[k]] .= fdx #ok
# 	end
# 	H_[solver.νinds_p[N-1], solver.xinds[N]] .= -Diagonal(@SVector ones(n)) #ok
#
# 	for k = 1:N-1
# 		ix = Z[k]._x #ok
# 		iu = Z[k]._u #ok
# 		fdu = ∇F[k][ix,iu] #ok
# 		H_[solver.νinds_p[k], solver.uinds[k]] .= fdu #okkk
# 	end
# 	return nothing
# end

function update_H_!(solver::DirectGamesSolver)
	n,m,N = size(solver)
	In = SMatrix{n,n}(Diagonal(SVector{n}(ones(n))))
	solver.H_ .*= 0.0#########################################################################################
	update_H_!(solver, In)
	return nothing
end

function update_H_!(solver::DirectGamesSolver, In::SMatrix{ln,ln,T,lnn}) where{ln,lnn,T}
	n,m,N = size(solver)
	n,m,pu,p = size(solver.model)

	Z = solver.Z
	H_ = solver.H_
	∇F = solver.∇F
	C = solver.C

	# Fill H_
	# Upper-left block
	for i = 1:p
		H_[solver.uinds_p[i][1], solver.uinds_H[i][1]] .= C[i].uu[1] #ok
		for k = 2:N-1
			H_[solver.xinds_p[i][k], solver.xinds[k]] .= C[i].xx[k] #ok
			H_[solver.uinds_p[i][k], solver.xinds[k]] .= C[i].ux[k] #ok
			H_[solver.xinds_p[i][k], solver.uinds_H[i][k]] .= C[i].ux[k]' #ok
			H_[solver.uinds_p[i][k], solver.uinds_H[i][k]] .= C[i].uu[k] #ok
		end
		H_[solver.xinds_p[i][N], solver.xinds[N]] .= C[i].xx[N] #ok
	end

			# # Dynamics Constraints
			# for k = 1:N-1
			# 	ix = Z[k]._x
			# 	iu = Z[k]._u
			# 	fdx = ∇F[k][ix,ix] # Ak [n,n]
			# 	fdu = ∇F[k][ix,iu] # Bk [n,m]
			# 	fdui = ∇F[k][ix,iu][:,pu[i]] # Bki [n,mi]
			# 	Iγi = Diagonal(solver.γ[i][k])
			# 	if k >= 2
			# 		H_[solver.xinds_p[i][k]  , solver.xinds[k]  ] .+=  fdx' *Iγi*fdx
			# 		H_[solver.xinds_p[i][k]  , solver.uinds[k]  ] .+=  fdx' *Iγi*fdu
			# 		H_[solver.xinds_p[i][k]  , solver.xinds[k+1]] .+= -fdx' *Iγi
			#
			# 		H_[solver.uinds_p[i][k]  , solver.xinds[k]  ] .+=  fdui'*Iγi*fdx
			#
			# 		H_[solver.xinds_p[i][k+1], solver.xinds[k]  ] .+= -      Iγi*fdx
			# 	end
			# 	H_[solver.uinds_p[i][k]  , solver.uinds[k]  ] .+=  fdui'*Iγi*fdu
			# 	H_[solver.uinds_p[i][k]  , solver.xinds[k+1]] .+= -fdui'*Iγi
			#
			# 	H_[solver.xinds_p[i][k+1], solver.uinds[k]  ] .+= -      Iγi*fdu
			# 	H_[solver.xinds_p[i][k+1], solver.xinds[k+1]] .+=        Iγi
			# end


			# for con in solver.dyn_constraints.constraints
			# 	for k in con.inds
			# 		rel_zind_i = rel_zinds(con,k,i,n,m,pu,p,N)
			# 		rel_zind = rel_zinds(con,k,n,m,pu,p,N)
			# 		zind_i = zinds(solver,con,k,i)
			# 		zind = zinds(solver,con,k)
			# 		∇c = con.∇c[k]
			# 		∇ci = con.∇c[k][:,rel_zind]
			# 		Iγi = Diagonal(solver.γ[i][k])
			# 		H_[zind_i,zind] += ∇ci'*Iγi*∇c# need check
	        #     end
	        # end ######## need to check all indices
		#
		# end
	for c in eachindex(solver.constraints.constraints)
		con = solver.constraints.constraints[c]
	# for con in solver.constraints.constraints
		update_H_con!(solver,con,n,m,p,N)
	end

	# Upper-right block
	for i = 1:p
		spui = solver.spu[i]
		for k = 2:N-1
			ix = Z[k]._x #ok
			fdx = ∇F[k][ix,ix] #ok
			H_[solver.xinds_p[i][k], solver.νinds[i][k]] .= fdx' #ok
			H_[solver.xinds_p[i][k], solver.νinds[i][k-1]] .= -In #ok
		end
		H_[solver.xinds_p[i][N], solver.νinds[i][N-1]] .= -In #ok

		for k = 1:N-1
			ix = Z[k]._x #ok
			iu = Z[k]._u #ok
			fdu = ∇F[k][ix,iu][:,spui] #ok
			H_[solver.uinds_p[i][k], solver.νinds[i][k]] .= fdu' #ok
		end
	end

	# Lower-left block
	for k = 2:N-1
		ix = Z[k]._x #ok
		fdx = ∇F[k][ix,ix] #ok
		H_[solver.νinds_p[k-1], solver.xinds[k]] .= -In #ok
		H_[solver.νinds_p[k], solver.xinds[k]] .= fdx #ok
	end
	H_[solver.νinds_p[N-1], solver.xinds[N]] .= -In #ok

	for k = 1:N-1
		ix = Z[k]._x #ok
		iu = Z[k]._u #ok
		fdu = ∇F[k][ix,iu] #ok
		H_[solver.νinds_p[k], solver.uinds[k]] .= fdu #okkk
	end
	return nothing
end

function update_H_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Stage},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	xinds = solver.xinds
	uinds = solver.uinds
	xinds_p = solver.xinds_p
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k == 1
				rel_zind_i = n .+ spu[i]
				rel_zind = n .+ sinds.sm
				zind_i = uinds_p[i][k]
				zind = uinds[k]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
				# set_view_H!(solver.stage_view_H[i][k],∇ci',Iμ,∇c)
			elseif k == N
				rel_zind_i = sinds.sn
				rel_zind = sinds.sn
				zind_i = xinds_p[i][k]
				zind = xinds[k]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
				# set_view_H!(solver.stage_view_H[i][k],∇ci',Iμ,∇c)
			else
				rel_zind_i = [sinds.sn; n.+ spu[i]]
				rel_zind = sinds.snm
				zind_i = [xinds_p[i][k]; uinds_p[i][k]]
				zind = [xinds[k]; uinds[k]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]

				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
				# set_mat_H!(solver.H_[zind_i,zind],∇ci',Iμ,∇c)
				# set_view_H!(solver.stage_view_H[i][k-1], ∇ci', Iμ, ∇c)
			end
		end
	end
	return nothing
end

function update_H_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State},
	n::Int, m::Int, p::Int, N::Int) where {T}
	xinds_p = solver.xinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k >= 2
				rel_zind_i = sinds.sn
				rel_zind = sinds.sn
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind] #ok
				Iμ = TO.penalty_matrix(con,j)
				set_view_H!(solver.state_view_H[i][k-1], ∇ci', Iμ, ∇c)
			end
		end
	end
	return nothing
end

function update_H_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Control},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	uinds = solver.uinds
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k <= N-1
				rel_zind_i = spu[i]
				rel_zind = sinds.sm
				# zind_i = uinds_p[i][k]
				# zind = uinds[k]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
				set_view_H!(solver.control_view_H[i][k], ∇ci', Iμ, ∇c)
			end
		end
	end
	return nothing
end

function update_H_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Coupled},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	xinds = solver.xinds
	uinds = solver.uinds
	xinds_p = solver.xinds_p
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k == 1
				rel_zind_i = [n .+ spu[i]; (n+m) .+ sinds.sn; (2*n+m) .+ spu[i]]
				rel_zind = [n .+ sinds.sm; (n+m) .+ sinds.snm]
				zind_i = [uinds_p[i][k]; xinds_p[i][k+1]; uinds_p[i][k+1]]
				zind = [uinds[k]; xinds[k+1]; uinds[k+1]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
			elseif k == N-1
				rel_zind_i = [sinds.sn; n .+ spu[i]; (n+m) .+ sinds.sn]
				rel_zind = sinds.s2nm
				zind_i = [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]]
				zind = [xinds[k]; uinds[k]; xinds[k+1]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
			elseif k == N
				rel_zind_i = sinds.sn
				rel_zind = sinds.sn
				zind_i = [xinds_p[i][k]]
				zind = [xinds[k]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
			else
				rel_zind_i = [sinds.sn; n .+ spu[i]; (n+m) .+ sinds.sn; (2*n+m) .+ spu[i]]
				rel_zind = sinds.s2n2m
				# zind_i = [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]; uinds_p[i][k+1]]
				# zind = [xinds[k]; uinds[k]; xinds[k+1]; uinds[k+1]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				# solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
				set_view_H!(solver.coupled_view_H[i][k-1], ∇ci', Iμ, ∇c)
			end
		end
	end
	return nothing
end

function update_H_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Dynamical},
	n::Int, m::Int, p::Int, N::Int) where {T}
	spu = solver.spu
	xinds = solver.xinds
	uinds = solver.uinds
	xinds_p = solver.xinds_p
	uinds_p = solver.uinds_p
	sinds = solver.sinds
	for i=1:p
		for (j,k) in enumerate(con.inds)
			if k == 1
				rel_zind_i = [n .+ spu[i]; (n+m) .+ sinds.sn]
				rel_zind = n .+ sinds.snm
				zind_i = [uinds_p[i][k]; xinds_p[i][k+1]]
				zind = [uinds[k]; xinds[k+1]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
			elseif k == N
				rel_zind_i = sinds.sn
				rel_zind = sinds.sn
				zind_i = [xinds_p[i][k]]
				zind = [xinds[k]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
			else
				rel_zind_i = [sinds.sn; n .+ spu[i]; (n+m) .+ sinds.sn]
				rel_zind = sinds.s2nm
				# zind_i = [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]]
				# zind = [xinds[k]; uinds[k]; xinds[k+1]]
				∇ci = con.∇c[j][:,rel_zind_i]
				∇c = con.∇c[j][:,rel_zind]
				Iμ = TO.penalty_matrix(con,j)
				# solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
				set_view_H!(solver.dynamical_view_H[i][k-1], ∇ci', Iμ, ∇c)
			end
		end
	end
	return nothing
end

function set_view_H!(v::SubArray, m1::SMatrix, m2::Diagonal{T,SVector{c,T}}, m3::SMatrix) where{c,T}
	v .+= m1*m2*m3
	return nothing
end
