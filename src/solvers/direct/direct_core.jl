export
	dual_ascent!,
	penalty_update!,
	inner_step!,
	line_search!,
	primal_dual_copy_update!,
	primal_dual_update!,
	indiv_cost

function dual_ascent!(solver::DirectGamesSolver)
	# @show "dual ascent"
    conSet = solver.constraints
	@show "here 5"
    for i in eachindex(conSet.constraints)
        TO.dual_update!(conSet.constraints[i])
    end
	@show "rrrr 5"
end

function penalty_update!(solver::DirectGamesSolver)
	# @show "penalty update"
	@show "here 6"
    for i in eachindex(solver.constraints.constraints)
        TO.penalty_update!(solver.constraints.constraints[i])
    end
	@show "rrrr 6"
	@show "here 7"
	for i in eachindex(solver.dyn_constraints.constraints)
		TO.penalty_update!(solver.dyn_constraints.constraints[i])
	end
	@show "rrrr 7"
end

function inner_step!(solver::DirectGamesSolver)
	# @show "inner step"
	solver.δY .= -(lu(solver.H_)\solver.g_) # better
	# solver.δY .= -solver.H_\solver.g_
	α, optimality_merit = line_search!(solver)
	record_inner_iteration!(solver, optimality_merit, α)
	primal_dual_update!(solver,α)
	dual_ascent!(solver)
	return nothing
end

function line_search!(solver::DirectGamesSolver)
	# @show "line search"
	α = solver.opts.α_init
	β = solver.opts.β
	τ = solver.opts.τ
	merit = mean(abs.(solver.g_))
	new_merit = 0.0
	for j = 1:solver.opts.iterations_linesearch
		primal_dual_copy_update!(solver, α)
		cost_expansion(solver.C, solver.obj, solver.Z̄, solver.model.pu, solver.model.p)
	    TO.evaluate!(solver.constraints, solver.Z̄)
	    TO.jacobian!(solver.constraints, solver.Z̄)
		TO.update_active_set!(solver.constraints,solver.Z̄)
		TO.evaluate!(solver.dyn_constraints, solver.Z̄)
		TO.discrete_jacobian!(solver.∇F, solver.model, solver.Z̄)

		update_g_!(solver, use_copy=true)
		new_merit = mean(abs.(solver.g_))
		if new_merit > (1-α*β)*merit
			α *= τ
		else
			regularization_update!(solver, :decrease)
			break
		end
		if j == solver.opts.iterations_linesearch
			regularization_update!(solver, :increase) ##############
			# solver.η[1] += solver.opts.bp_reg_fp # Very harmful to convergence
			solver.stats.dJ_zero_counter += 1
		end
	end
	return α, new_merit
end

function primal_dual_copy_update!(solver::DirectGamesSolver, α::T) where {T}
	# # @show "primal_dual copy update"
	n,m,N = size(solver)
	n,m,pu,p = size(solver.model)
	# handles u1
	solver.Z̄[1].z = solver.Z[1].z + α * [@SVector zeros(n); solver.δY[solver.uinds[1]]]
	for k = 2:N-1
		solver.Z̄[k].z = solver.Z[k].z + α * solver.δY[[solver.xinds[k]; solver.uinds[k]]]
	end
	# handles xN
	solver.Z̄[N].z = solver.Z[N].z + α * [solver.δY[solver.xinds[N]]; @SVector zeros(m)]

	for i = 1:p
		for k = 1:N-1
			solver.ν_[i][k] = solver.ν[i][k] + α * solver.δY[solver.νinds[i][k]]
		end
	end
	return nothing
end

function primal_dual_update!(solver::DirectGamesSolver, α::T) where {T}
	# @show "primal dual update"
	n,m,N = size(solver)
	n,m,pu,p = size(solver.model)
	# handles u1
	solver.Z[1].z += α * [@SVector zeros(n); solver.δY[solver.uinds[1]]]
	for k = 2:N-1
		solver.Z[k].z += α * solver.δY[[solver.xinds[k]; solver.uinds[k]]]
	end
	# handles xN
	solver.Z[N].z += α * [solver.δY[solver.xinds[N]]; @SVector zeros(m)]

	for i = 1:p
		for k = 1:N-1
			solver.ν[i][k] += α * solver.δY[solver.νinds[i][k]]
		end
	end
	return nothing
end

function indiv_cost(solver::DirectGamesSolver)
	# @show "indiv cost"
	_J = TO.get_J.(solver.obj)
	J = sum.(_J)
	return J
end












#
# function set_state_view_H!(
# 	state_view::SubArray,
# 	m1::SMatrix,
# 	m2,
# 	m3::SMatrix,
# 	)
# 	state_view .+= m1*m2*m3
# 	return nothing
# end
#
# function set_mat_H!(
# 	mat::SparseMatrixCSC,
# 	m1::SMatrix,
# 	m2,
# 	m3::SMatrix,
# 	)
# 	mat += m1*m2*m3
# 	return nothing
# end


# function update_g_con_old!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State},
# 	n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T}
# 	for i=1:p
# 		for k in intersect(con.inds,2:N)
# 			# rel_zind_i = rel_zinds(con,solver.sinds_p,k,i,N)
# 			# @show rel_zind_i
# 			zind_i = zinds(solver,con,k,i)
# 			# ∇ci = con.∇c[k][:,rel_zind_i]
# 			∇ci = con.∇c[k][:,[solver.Z[k]._x; (-6 .+ solver.Z[k]._u)]]
# 			# ∇ci = con.∇c[k][:,s]
# 			Iμ = TO.penalty_matrix(con,k)
# 			# @inbounds solver.g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] # best
# 		end
# 	end
# 	return nothing
# end
#
# function update_g_con_old!(solver::DirectGamesSolver{T}, con::ConstraintVals{T},
# 	n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T}
# 	# println("Not state constraint")
# 	for i=1:p
# 		for k in con.inds
# 			rel_zind_i = rel_zinds(con,solver.sinds_p,k,i,N)
# 			zind_i = zinds(solver,con,k,i)
# 			∇ci = con.∇c[k][:,rel_zind_i]
# 			Iμ = TO.penalty_matrix(con,k)
# 			@inbounds solver.g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k]
# 		end
# 	end
# 	return nothing
# end
# a = 10



# function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State},
# 	n::Int, m::Int, p::Int, N::Int) where {T}
# 	sinds = solver.sinds
# 	for i=1:p
# 		for (j,k) in enumerate(con.inds)
# 			if k >= 2
# 				rel_zind_i = sinds.sn
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			end
# 		end
# 	end
# 	return nothing
# end
#
# function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Control},
# 	n::Int, m::Int, p::Int, N::Int) where {T}
# 	sinds = solver.sinds
# 	for i=1:p
# 		for (j,k) in enumerate(con.inds)
# 			if k <= N-1
# 				rel_zind_i = sinds.sm
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			end
# 		end
# 	end
# 	return nothing
# end
#
# function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Coupled},
# 	n::Int, m::Int, p::Int, N::Int) where {T}
# 	sinds = solver.sinds
# 	for i=1:p
# 		for (j,k) in enumerate(con.inds)
# 			if k == 1
# 				rel_zind_i = [n .+ sinds.sm; (n+m) .+ sinds.snm]
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			elseif k == N-1
# 				rel_zind_i = sinds.s2nm
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			elseif k == N
# 				rel_zind_i = sinds.sn
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			else
# 				rel_zind_i = sinds.s2n2m
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			end
# 		end
# 	end
# 	return nothing
# end
#
# function update_g_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Dynamical},
# 	n::Int, m::Int, p::Int, N::Int) where {T}
# 	sinds = solver.sinds
# 	for i=1:p
# 		for (j,k) in enumerate(con.inds)
# 			if k == 1
# 				rel_zind_i = n .+ sinds.snm
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			elseif k == N
# 				rel_zind_i = sinds.sn
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			else
# 				rel_zind_i = sinds.s2nm
# 				∇ci = con.∇c[j][:,rel_zind_i]
# 				zind_i = zinds(solver,con,k,i)
# 				Iμ = TO.penalty_matrix(con,j)
# 				@inbounds solver.g_[zind_i] += ∇ci'*con.λ[j] + ∇ci'*Iμ*con.vals[j] # best
# 			end
# 		end
# 	end
# 	return nothing
# end

# k = N
# i = 2
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames,k,i,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames,k,i,n,m,N)
# @allocated rel_zind_i = con_rel_zinds(con,solver_directgames,k,i,n,m,N)
# test_allocation(rel_zinds, (con,solver_directgames,k,i,N))

# p = 3
# N = 41
# con = solver_directgames.constraints.constraints[1]
# # con = solver_directgames.constraints.constraints[end]
# # for con in solver_directgames.constraints.constraints
# # 	println(typeof(con))
# # end
# # @btime update_g_con!($solver_directgames, $con, $n,$m,$pu,$p, $N)
# # @btime update_g_con!($solver_directgames, $con, $n,$m,$pu,$p, $N)
# @allocated update_g_con!(solver_directgames, con, n,m,pu,p, N)
# @allocated update_g_con!(solver_directgames, con, n,m,pu,p, N)
# @allocated update_g_con!(solver_directgames, con, n,m,pu,p, N)
# @allocated update_g_con!(solver_directgames, con, n,m,pu,p, N)
# test_allocation(update_g_con!,(solver_directgames, con, n,m,pu,p, N))

# a = 10
# a = 10
# a = 10
# a = 10
# a = 10
#


#
# @allocated update_g_!(solver_directgames)
# @allocated update_g_!(solver_directgames)
# @allocated update_g_!(solver_directgames)
# @allocated update_g_!(solver_directgames)
# @allocated update_g_!(solver_directgames)
# @btime update_g_!(solver_directgames)
# @btime update_g_old!(solver_directgames)
# #
#
#

# con = solver_directgames.constraints.constraints[1]
# # con = solver_directgames.constraints.constraints[end]
# # for con in solver_directgames.constraints.constraints
# # 	println(typeof(con))
# # end
# # @btime update_g_con!($solver_directgames, $con, $n,$m,$pu,$p, $N)
# # @btime update_g_con!($solver_directgames, $con, $n,$m,$pu,$p, $N)
# @allocated update_g_con!(solver_directgames, con, n,m,p, N)
# @allocated update_g_con!(solver_directgames, con, n,m,p, N)
# @allocated update_g_con!(solver_directgames, con, n,m,p, N)
# @allocated update_g_con!(solver_directgames, con, n,m,p, N)
# test_allocation(update_g_con!,(solver_directgames, con, n,m,p, N))


#
# con = solver_directgames.constraints.constraints[1]
# con = solver_directgames.constraints.constraints[end]
# @allocated update_H_con!(solver_directgames, con, n,m,p, N)
# @allocated update_H_con!(solver_directgames, con, n,m,p, N)
# @allocated update_H_con!(solver_directgames, con, n,m,p, N)
# @allocated update_H_con!(solver_directgames, con, n,m,p, N)
# # @btime update_H_con!(solver_directgames, con, n,m,p, N)
# test_allocation(update_H_con!,(solver_directgames, con, n,m,p, N))
#
#
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @btime update_H_!(solver_directgames)
# @allocated update_H_old!(solver_directgames)
# @btime update_H_old!(solver_directgames)
# #
# #

# function set_state_view_H2!(
# 	state_view::SubArray,
# 	m1::SMatrix,
# 	m2,
# 	m3::SMatrix,
# 	)
# 	state_view .+= m1*m2*m3
# 	return nothing
# end

# sm1 = SMatrix{n,n}(ones(n,n))
# sm2 = Diagonal(SVector{n}(ones(n)))
# sm3 = SMatrix{n,n}(ones(n,n))
# @allocated set_state_view_H2!(solver_directgames.state_view_H[1][1], sm1,sm2,sm3)
# @allocated set_state_view_H2!(solver_directgames.state_view_H[1][1], sm1,sm2,sm3)
# @allocated set_state_view_H2!(solver_directgames.state_view_H[1][1], sm1,sm2,sm3)
# @allocated set_state_view_H2!(solver_directgames.state_view_H[1][1], sm1,sm2,sm3)


# sm1 = SMatrix{n,n}(ones(n,n))
# sm2 = Diagonal(SVector{n}(ones(n)))
# sm3 = SMatrix{n,n}(ones(n,n))
# viewH = view(solver_directgames.H_, 1:12, 1:12)
# @allocated set_mat_H!(viewH, sm1,sm2,sm3)
# @allocated set_mat_H!(viewH, sm1,sm2,sm3)
# @allocated set_mat_H!(viewH, sm1,sm2,sm3)
# @allocated set_mat_H!(solver_directgames.H_[1:12, 1:12], sm1,sm2,sm3)
# @allocated set_mat_H!(solver_directgames.H_[1:12, 1:12], sm1,sm2,sm3)
# @allocated set_mat_H!(solver_directgames.H_[1:12, 1:12], sm1,sm2,sm3)
#

#
#
# function update_H_con!(solver::DirectGamesSolver, con::ConstraintVals{T,State},
# 	n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T}
# 	for i=1:p
# 		for k in intersect(con.inds,2:N)
# 			rel_zind_i = rel_zinds(con,solver.sinds,k,i,N)
# 			rel_zind = rel_zinds(con,solver.sinds,k,N)
# 			∇ci = con.∇c[k][:,rel_zind_i] #ok
# 			∇c = con.∇c[k][:,rel_zind] #ok
# 			Iμ = TO.penalty_matrix(con,k)
# 			set_state_view_H!(solver.state_view_H[i][k-1], ∇ci', Iμ, ∇c)
# 		end
# 	end
# 	return nothing
# end
#
# function update_H_con!(solver::DirectGamesSolver{T}, con::ConstraintVals{T},
# 	n::Int, m::Int, pu::Vector{Vector{Int}}, p::Int, N::Int) where {T}
# 	# println("Not state constraint")
# 	for i=1:p
# 		for k in con.inds
# 			rel_zind_i = rel_zinds(con,solver.sinds_p,k,i,N)
# 			rel_zind = rel_zinds(con,solver.sinds,k,N)
# 			zind_i = zinds(solver,con,k,i)
# 			zind = zinds(solver,con,k)
# 			∇ci = con.∇c[k][:,rel_zind_i] #ok
# 			∇c = con.∇c[k][:,rel_zind] #ok
# 			Iμ = TO.penalty_matrix(con,k)
# 			solver.H_[zind_i,zind] .+= ∇ci'*Iμ*∇c
# 		end
# 	end
# 	return nothing
# end

#
#
# @allocated update_H_!(solver_directgames, Im)
# @allocated update_H_!(solver_directgames, Im)
# @allocated update_H_!(solver_directgames, Im)
# @allocated update_H_!(solver_directgames, Im)
# @allocated update_H_!(solver_directgames, Im)
# @btime update_H_!(solver_directgames, Im)
# @btime update_H_old!(solver_directgames)


# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @allocated update_H_!(solver_directgames)
# @btime update_H_!(solver_directgames)
# @btime update_H_old!(solver_directgames)


#
#

# p = 3
# N = 41
# con = solver_directgames.constraints.constraints[1]
# @allocated update_H_con!(solver_directgames, con, n,m,pu,p, N)
# @allocated update_H_con!(solver_directgames, con, n,m,pu,p, N)
# @allocated update_H_con!(solver_directgames, con, n,m,pu,p, N)
# @allocated update_H_con!(solver_directgames, con, n,m,pu,p, N)
# test_allocation(update_H_con!,(solver_directgames, con, n,m,pu,p, N))


#
# @allocated update_H_con!(solver_directgames, con, n,m,pu,p, N,temp1, temp2)
# @btime update_H_con!(solver_directgames, con, n,m,pu,p, N,temp1, temp2)

#
# function set_state_view_H!(
# 	state_view::SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{n,Int},SVector{n,Int}},false},
# 	m1::SMatrix{n,c,T,nc},
# 	m2::SMatrix{c,c,T,cc},
# 	m3::SMatrix{c,n,T,nc},
# 	) where {n,c,nc,cc,T}
# 	state_view .+= m1*m2*m3
# 	return nothing
# end


#
#
# mm1 = SMatrix{12,1}(ones(12,1))
# mm2 = SMatrix{1,1}(Diagonal(SVector{1}(ones(1))))
# mm3 = SMatrix{1,12}(ones(1,12))
# # set_state_view_H!(solver_directgames.state_view_H[1][1], mm1,mm2,mm3)
# @allocated set_state_view_H_clean!(solver_directgames.state_view_H[1][1], mm1,mm2,mm3)
# @allocated set_state_view_H_clean!(solver_directgames.state_view_H[1][1], mm1,mm2,mm3)
# @allocated set_state_view_H_clean!(solver_directgames.state_view_H[1][1], mm1,mm2,mm3)
# @allocated set_state_view_H_clean!(solver_directgames.state_view_H[1][1], mm1,mm2,mm3)
#
#
# function set_view!(spMat::SparseMatrixCSC{T,Int}, mat::SMatrix{ln,lm,T,ll}) where {ln,lm,ll,T}
# 	spMat .= mat
# 	return nothing
# end
#
# function h_dumb5(
# 	H::SparseMatrixCSC{T,Int},
# 	# sview::SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{n,Int},SVector{n,Int}},false},
# 	temp1::SMatrix{ln,lm,T,ll},
# 	temp2::SMatrix{ln,lm,T,ll},
# 	temp3::SMatrix{ln,lm,T,ll},
# 	) where {ln,lm,ll
# 	,T}
# 	# H[1:12, 1:12] .= 0.0
# 	for i = 1:100
# 		H .+= temp1*temp2*temp3
# 		# sview .+= temp1*temp2*temp3
# 	end
# 	return nothing
# end
#
# temp = SMatrix{12,12}(ones(12,12))
# hh = solver_directgames.H_[1:12,1:12]
# @allocated h_dumb5(hh, temp, temp, temp)
# @allocated h_dumb5(hh, temp, temp, temp)
# @allocated h_dumb5(hh, temp, temp, temp)
# @allocated h_dumb5(hh, temp, temp, temp)
# @btime h_dumb5(hh, temp, temp, temp)
#
# function h_dumb4(
# 	# H::SparseMatrixCSC{T,Int},
# 	sview::SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{n,Int},SVector{n,Int}},false},
# 	temp1::SMatrix{ln,lm,T,ll},
# 	temp2::SMatrix{ln,lm,T,ll},
# 	temp3::SMatrix{ln,lm,T,ll},
# 	) where {ln,lm,ll
# 	,T}
# 	# H[1:12, 1:12] .= 0.0
# 	for i = 1:100
# 		# H[1:12, 1:12] .+= temp1*temp2
# 		sview .+= temp1*temp2*temp3
# 	end
# 	return nothing
# end
#
# temp = SMatrix{12,12}(ones(12,12))
# hh = solver_directgames.H_[1:12,1:12]
# hh = solver_directgames.state_view_H[1][1]
# @allocated h_dumb4(hh, temp, temp, temp)
# @allocated h_dumb4(hh, temp, temp, temp)
# @allocated h_dumb4(hh, temp, temp, temp)
# @allocated h_dumb4(hh, temp, temp, temp)
# @btime h_dumb4(hh, temp, temp, temp)
#
#





#
#
#
# Im = SMatrix{n,n}(Diagonal(SVector{n}(ones(n))))
#
#
#
# function h_dumb(H::SparseMatrixCSC{T,Int}, In::Diagonal{Float64,SVector{nn,T}}) where {nn,T}
# 	# H[1:12, 1:12] .= 0.0
# 	for i = 1:100
# 		H[1:12, 1:12] .= In
# 	end
# 	return nothing
# end
#
# Inn = Diagonal(SVector{n}(ones(n)))
# @allocated h_dumb(solver_directgames.H_, Inn)
# @allocated h_dumb(solver_directgames.H_, Inn)
# @allocated h_dumb(solver_directgames.H_, Inn)
# @allocated h_dumb(solver_directgames.H_, Inn)
# @btime h_dumb(solver_directgames.H_, Inn)
#
#
# function h_dumb2(H::SparseMatrixCSC{T,Int}, In::SVector{nn,T}) where {nn,T}
# 	inds = SVector{12}(diagind(H[1:12, 1:12]))
# 	for i=1:100
# 		H[1:12, 1:12] .= 0.0
# 		H[inds] .= In
# 	end
# 	return nothing
# end
#
# Iv = SVector{n}(ones(n))
# @allocated h_dumb2(solver_directgames.H_, Iv)
# @allocated h_dumb2(solver_directgames.H_, Iv)
# @allocated h_dumb2(solver_directgames.H_, Iv)
# @allocated h_dumb2(solver_directgames.H_, Iv)
# @btime h_dumb2(solver_directgames.H_, Iv)
#
#
#
# function h_dumb3(H::SparseMatrixCSC{T,Int}, In::SMatrix{nn,nn,T,NN}) where {nn,NN,T}
# 	# Imm = SMatrix{n,n}(Matrix{T}(I,n,n))
# 	for i=1:100
# 		H[1:12, 1:12] .= 0.0
# 		for j=1:12
# 			H[j,j] = 1.0
# 		end
# 	end
# 	return nothing
# end
#
# Im = SMatrix{n,n}(ones(n,n))
# @allocated h_dumb3(solver_directgames.H_, Im)
# @allocated h_dumb3(solver_directgames.H_, Im)
# @allocated h_dumb3(solver_directgames.H_, Im)
# @allocated h_dumb3(solver_directgames.H_, Im)
# @btime h_dumb3(solver_directgames.H_, Im)
#
#
#
#
#
#
# A = ones(12,12)
# inds = [[i,i] for i =1:12]
# A[diagind(A)] = 10*ones(12)
# A
# diagind(A)
#


#
#
#
# ############################################
# #############################################
# #############################################
# ##############################################
# ###########################################
# ###############################################
# ######################
#



#
# solver_directgames.g_ .*= 0.0
# maximum(abs.(solver_directgames.g_))
# update_g_!(solver_directgames)
# maximum(abs.(solver_directgames.g_))
# new = copy(solver_directgames.g_)
#
# solver_directgames.g_ .*= 0.0
# maximum(abs.(solver_directgames.g_))
# update_g_old!(solver_directgames)
# maximum(abs.(solver_directgames.g_))
# old = copy(solver_directgames.g_)
#
#
# maximum(abs.(old - new))
# a = 10
#

# #
# # function update_g_lite!(solver::DirectGamesSolver; use_copy=false::Bool)
# # 	# println("************************GGGGGGGGGGGGGGGGGGGGGGGGGGG")
# # 	n,m,N = size(solver)
# # 	n,m,pu,p = size(solver.model)
# # 	if use_copy
# # 		Z = solver.Z̄
# # 		ν = solver.ν_
# # 	else
# # 		Z = solver.Z
# # 		ν = solver.ν
# # 	end
# #
# # 	g_ = solver.g_
# # 	∇F = solver.∇F
# # 	C = solver.C
# # 	dyn = solver.dyn_constraints.constraints[1].vals
# # 	# dyn = Vector{SVector{n,T}}(solver.dyn_constraints.constraints[1].vals)
# #
# #
# #
# # 	# # spu = [SVector{length(pu[i])}(pu[i]) for i in 1:length(pu)]
# # 	# # spu = [SVector{2}([1,2]),SVector{2}([3,4]),SVector{2}([5,6])]
# #
# # 	#
# # 	# # Fill g_
# # 	# for i = 1:p::Int
# # 	# 	# # Dynamics and cost
# # 	# 	# # uk
# # 	# 	# spui = solver.spu[i]
# # 	# 	# for k = 1:N-1
# # 	# 	# 	ix = Z[k]._x
# # 	# 	# 	iu = Z[k]._u
# # 	# 	# 	fdu = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
# # 	# 	# 	g_[solver.uinds_p[i][k]] = C[i].u[k] + fdu'*ν[i][k] #ok
# # 	# 	# end
# # 	# 	# # xk
# # 	# 	# for k = 2:N-1
# # 	# 	# 	ix = Z[k]._x
# # 	# 	# 	fdx = ∇F[k][ix,ix] # Ak [n,n]
# # 	# 	# 	g_[solver.xinds_p[i][k]] = C[i].x[k] + fdx'*ν[i][k] - ν[i][k-1] #ok
# # 	# 	# end
# # 	# 	# g_[solver.xinds_p[i][N]] = C[i].x[N] - ν[i][N-1] #ok
# # 	#
# # 	# 	# # Constraints
# # 	# 	# for con in solver.constraints.constraints
# # 	# 	# 	if typeof(con).parameters[2] == State
# # 	# 	# 		for k in intersect(con.inds,2:N)
# # 	# 	# 			# rel_zind_i = rel_zinds(con,solver.Z[k],k,i,n,m,pu,p,N)
# # 	# 	# 			# zind_i = zinds(solver,con,k,i)
# # 	# 	# 			# ∇ci = con.∇c[k][:,rel_zind_i]
# # 	# 	# 			# Iμ = TO.penalty_matrix(con,k)
# # 	# 	# 			# @inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] # best
# # 	# 	# 			# view_g = solver.state_view_g[i][k-1] #bad
# # 	# 	# 			# view_g .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #bad
# # 	# 	# 			# solver.state_view_g[i][k-1] .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #Bad
# # 	# 	# 		end
# # 	# 	# 	else
# # 	# 	# 		for k in con.inds
# # 	# 	# 			# rel_zind_i = rel_zinds(con,solver.Z[k],k,i,n,m,pu,p,N)
# # 	# 	# 			# zind_i = zinds(solver,con,k,i)
# # 	# 	# 			# ∇ci = con.∇c[k][:,rel_zind_i]
# # 	# 	# 			# Iμ = TO.penalty_matrix(con,k)
# # 	# 	# 			# @inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k]
# # 	# 	# 		end
# # 	# 	# 	end
# # 	# 	# end
# # 	#
# # 	# 	# Dynamics constraints
# # 	# 	# for k = 1:N-1
# # 	# 	# 	Iγi = Diagonal(solver.γ[i][k])
# # 	# 	# 	ix = Z[k]._x
# # 	# 	# 	iu = Z[k]._u
# # 	# 	# 	fdui = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
# # 	# 	# 	# g_[solver.uinds_p[i][k]] += fdui'*Iγi*dyn[k]# need check
# # 	# 	# end
# # 	# 	# for k = 2:N-1 ##### MAYBE K = 3:N-1
# # 	# 	# 	Iγi = Diagonal(solver.γ[i][k])
# # 	# 	# 	ix = Z[k]._x
# # 	# 	# 	fdx = ∇F[k][ix,ix] # Ak [n,n]
# # 	# 	# 	g_[solver.xinds_p[i][k]] += fdx'*Iγi*dyn[k]# need check
# # 	# 	# end
# # 	# 	# for k = 1:N-1
# # 	# 	# 	Iγi = Diagonal(solver.γ[i][k])
# # 	# 	# 	g_[solver.xinds_p[i][k+1]] = dyn[k]# need check
# # 	# 	# 	# g_[solver.xinds_p[i][k+1]] = -Iγi*dyn[k]# need check
# # 	# 	# 	# g_[solver.xinds_p[i][k+1]] += -Iγi*dyn[k]# need check
# # 	# 	# end
# # 	# end
# #
# # 	# # Dynamics Constraints
# # 	# col_off = Int((n*p+m)*(N-1)) #ok
# # 	# d = dyn[1]
# # 	# for k = 1:N-1::Int
# # 	#
# # 	# 	ii = Z[k]._x
# # 	# 	# inds =  col_off+(k-1)*n .+ ix
# # 	# 	# @show typeof(inds)
# # 	# 	# g_[ii] = d #ok
# # 	# end
# #
# # 	# Dynamics Constraints
# # 	col_off = (n*p+m)*(N-1) #ok
# # 	# ii = solver_directgames.Z[1]._x
# # 	# inds = SVector{12}([i for i=1:12])
# # 	# ddd = SVector{12}([1.0 for i=1:12])
# # 	# ddd += dyn[1]
# # 	# @show typeof(SVector{12}([1.0 for i=1:12]))
# # 	# @show typeof(dyn[1])
# # 	# indss = solver.Z[1]._x
# #
# # 	# for k = 1:N-1::Int
# # 	# # 	# g_[col_off+(k-1)*n .+ (1:n)] = dyn[k] #ok
# # 	# # 	@show typeof(dyn[k])
# # 	# # 	@show typeof(g_[col_off+(k-1)*n .+ ii])
# # 	# # 	@show typeof(col_off+(k-1)*n .+ ii)
# # 	# # 	g_[col_off+(k-1)*n .+ ii] = dyn[k] #ok
# # 	# 	# g_[ii .+ (col_off+(k-1)*n)] = dyn[k]
# # 	# 	# g_[ii] = dyn[k]
# # 	# 	# g_[ii]
# # 	# 	indss = solver.Z[k]._x
# # 	# 	# ddd = dyn[k]
# # 	# 	g_[indss.+ (col_off+(k-1)*n)] += 100000*dyn[k]
# # 	# 	# g_[indss.+ (col_off+(k-1)*n)] .+= 1.0
# # 	# 	# g_[indss.+ (col_off+(k-1)*n)] = solver.dyn_constraints.constraints[1].vals[1]
# # 	# end
# #
# # 	# # # Dynamics Constraints
# # 	# col_off = (n*p+m)*(N-1) #ok
# # 	# for k = 1:N-1
# # 	# 	g_[col_off+(k-1)*n .+ (1:n)] = dyn[k] #ok
# # 	# end
# # 	update_g_lite_dyn!(solver, dyn, n, m, p, N)
# #
# # 	# g_dyn!(g_, n, N, indss, col_off, dyn)
# # 	# @show g_[12]
# # 	return nothing
# # end
# #
#
# # function update_g_lite_best!(solver::DirectGamesSolver; use_copy=false::Bool)
# # 	n,m,N = size(solver)
# # 	n,m,pu,p = size(solver.model)
# # 	if use_copy
# # 		Z = solver.Z̄
# # 		ν = solver.ν_
# # 	else
# # 		Z = solver.Z
# # 		ν = solver.ν
# # 	end
# #
# # 	g_ = solver.g_
# # 	∇F = solver.∇F
# # 	C = solver.C
# # 	dyn = solver.dyn_constraints.constraints[1].vals
# #
# # 	# Fill g_
# # 	for i = 1:p
# # 		# Dynamics and cost
# # 		# uk
# # 		spui = solver.spu[i]
# # 		for k = 1:N-1
# # 			ix = Z[k]._x
# # 			iu = Z[k]._u
# # 			fdu = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
# # 			g_[solver.uinds_p[i][k]] = C[i].u[k] + fdu'*ν[i][k] #ok
# # 		end
# # 		# xk
# # 		for k = 2:N-1
# # 			ix = Z[k]._x
# # 			fdx = ∇F[k][ix,ix] # Ak [n,n]
# # 			g_[solver.xinds_p[i][k]] = C[i].x[k] + fdx'*ν[i][k] - ν[i][k-1] #ok
# # 		end
# # 		g_[solver.xinds_p[i][N]] = C[i].x[N] - ν[i][N-1] #ok
# #
# # 		# # Constraints
# # 		# for con in solver.constraints.constraints
# # 		# 	if typeof(con).parameters[2] == State
# # 		# 		for k in intersect(con.inds,2:N)
# # 		# 			rel_zind_i = rel_zinds(con,solver.Z[k],k,i,n,m,pu,p,N)
# # 		# 			zind_i = zinds(solver,con,k,i)
# # 		# 			∇ci = con.∇c[k][:,rel_zind_i]
# # 		# 			Iμ = TO.penalty_matrix(con,k)
# # 		# 			@inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] # best
# # 		# 			# view_g = solver.state_view_g[i][k-1] #bad
# # 		# 			# view_g .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #bad
# # 		# 			# solver.state_view_g[i][k-1] .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #Bad
# # 		# 		end
# # 		# 	else
# # 		# 		for k in con.inds
# # 		# 			rel_zind_i = rel_zinds(con,solver.Z[k],k,i,n,m,pu,p,N)
# # 		# 			zind_i = zinds(solver,con,k,i)
# # 		# 			∇ci = con.∇c[k][:,rel_zind_i]
# # 		# 			Iμ = TO.penalty_matrix(con,k)
# # 		# 			@inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k]
# # 		# 		end
# # 		# 	end
# # 		# end
# # 		#
# # 		# # Dynamics constraints
# # 		# for k = 1:N-1
# # 		# 	Iγi = Diagonal(solver.γ[i][k])
# # 		# 	ix = Z[k]._x
# # 		# 	iu = Z[k]._u
# # 		# 	fdui = ∇F[k][ix,iu][:,pu[i]] # Bki [n,mi]
# # 		# 	g_[solver.uinds_p[i][k]] += fdui'*Iγi*dyn[k]# need check
# # 		# end
# # 		# for k = 2:N-1 ##### MAYBE K = 3:N-1
# # 		# 	Iγi = Diagonal(solver.γ[i][k])
# # 		# 	ix = Z[k]._x
# # 		# 	fdx = ∇F[k][ix,ix] # Ak [n,n]
# # 		# 	g_[solver.xinds_p[i][k]] += fdx'*Iγi*dyn[k]# need check
# # 		# end
# # 		for k = 1:N-1
# # 			Iγi = Diagonal(solver.γ[i][k])
# # 			# g_[solver.xinds_p[i][k+1]] += -Iγi*dyn[k]# need check
# # 		end
# # 	end
# # 	# # Dynamics Constraints
# # 	# col_off = (n*p+m)*(N-1) #ok
# # 	# for k = 1:N-1
# # 	# 	g_[col_off+(k-1)*n .+ (1:n)] = dyn[k] #ok
# # 	# end
# # 	update_g_lite_dyn!(solver, dyn, n, m, p, N)
# # 	return nothing
# # end
#
# function update_g_lite!(solver::DirectGamesSolver; use_copy=false::Bool)
# 	dyn = solver.dyn_constraints.constraints[1].vals
# 	update_g_lite!(solver,dyn)
# 	return nothing
# end
#
#
# function update_g_lite!(solver::DirectGamesSolver, dyn::Vector{SVector{nn,T}}; use_copy=false::Bool) where {nn,T}
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
#
# 	# Fill g_
# 	for con in solver.constraints.constraints
# 		update_g_con!(solver,con,n,m,pu,p,N)
# 	end
#
# 	for i = 1:p
# 		# Dynamics and cost
# 		# uk
# 		spui = solver.spu[i]
# 		for k = 1:N-1
# 			ix = Z[k]._x
# 			iu = Z[k]._u
# 			fdu = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
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
# 		# # Constraints
# 		# for con in solver.constraints.constraints
# 		# 	# @show typeof(con)
# 		# 	if typeof(con).parameters[2] == State
# 		# 		for k in intersect(con.inds,2:N)
# 		# 			# rel_zind_i = rel_zinds(con,solver.Z[k],k,i,n,m,pu,p,N)
# 		# 			# zind_i = zinds(solver,con,k,i)
# 		# 			# ∇ci = con.∇c[k][:,rel_zind_i]
# 		# 			# Iμ = TO.penalty_matrix(con,k)
# 		# 			# @inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] # best
# 		# 			# # view_g = solver.state_view_g[i][k-1] #bad
# 		# 			# # view_g .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #bad
# 		# 			# # solver.state_view_g[i][k-1] .+= ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k] #Bad
# 		# 		end
# 		# 	else
# 		# 		for k in con.inds
# 		# 			# rel_zind_i = rel_zinds(con,solver.Z[k],k,i,n,m,pu,p,N)
# 		# 		# 	zind_i = zinds(solver,con,k,i)
# 		# 		# 	∇ci = con.∇c[k][:,rel_zind_i]
# 		# 		# 	Iμ = TO.penalty_matrix(con,k)
# 		# 		# 	@inbounds g_[zind_i] += ∇ci'*con.λ[k] + ∇ci'*Iμ*con.vals[k]
# 		# 		end
# 		# 	end
# 		# end
#
# 		# Dynamics constraints
# 		for k = 1:N-1
# 			Iγi = Diagonal(solver.γ[i][k])
# 			ix = Z[k]._x
# 			iu = Z[k]._u
# 			fdui = ∇F[k][ix,iu][:,spui] # Bki [n,mi]
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
#
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @btime update_g_lite!($solver_directgames)
# @btime update_g_!($solver_directgames)
#
#
# p = 3
# N = 10
# con = solver_directgames.constraints.constraints[1]
# @allocated update_g_con!(solver_directgames, con, p, N)
# @allocated update_g_con!(solver_directgames, con, p, N)
# @allocated update_g_con!(solver_directgames, con, p, N)
# @allocated update_g_con!(solver_directgames, con, p, N)
# @btime update_g_con!(solver_directgames, con, p, N)
#
#
#
#
#
# dcon = solver_directgames.dyn_constraints.constraints[1]
# @allocated update_g_con!(solver_directgames, dcon, p, N)
# @allocated update_g_con!(solver_directgames, dcon, p, N)
# @allocated update_g_con!(solver_directgames, dcon, p, N)
# @allocated update_g_con!(solver_directgames, dcon, p, N)
# @btime update_g_con!(solver_directgames, dcon, p, N)
#
#
#
#
#
#
#
#
#
#
#
#
# function update_g_lite_dyn!(solver::DirectGamesSolver{T},dyn::Vector{SVector{nn,T}},n::Int,m::Int,p::Int,N::Int) where {nn,T}
# 	col_off = (n*p+m)*(N-1) #ok
# 	for k = 1:N-1
# 		solver.g_[col_off+(k-1)*n .+ (1:n)] = dyn[k] #ok
# 	end
# 	return nothing
# end
# @allocated update_g_lite_dyn!(solver_directgames,din,n,m,p,N)
# @allocated update_g_lite_dyn!(solver_directgames,din,n,m,p,N)
# @allocated update_g_lite_dyn!(solver_directgames,din,n,m,p,N)
# @allocated update_g_lite_dyn!(solver_directgames,din,n,m,p,N)
#
#
#
#
# # inds = SVector{12}([i for i=1:12])
#
# function g_dyn!(g_::Vector{T},
# 	nn::Int, N::Int, indss::SVector{n,Int}, col_off::Int, dyn::Vector{SVector{n,T}}) where {n,T}
# 	for k = 1:N-1::Int
# 		g_[indss.+ (col_off+(k-1)*nn)] = dyn[k]
# 	end
# 	return nothing
# end
# indss = solver_directgames.Z[1]._x
# col_off = 11
# @allocated g_dyn!(solver_directgames.g_, n, N, indss, col_off, din)
# @show solver_directgames.g_[1]
# @allocated g_dyn!(solver_directgames.g_, n, N, indss, col_off, din)
# @allocated g_dyn!(solver_directgames.g_, n, N, indss, col_off, din)
# @allocated g_dyn!(solver_directgames.g_, n, N, indss, col_off, din)
# @allocated g_dyn!(solver_directgames.g_, n, N, indss, col_off, din)
#
#
#
#
# struct Bundle2{T}
# 	g::Vector{T}
# 	din::Vector{SVector{12,T}}
# 	inds::SVector{12,Int}
# 	sv::SVector{12,T}
# 	N::Int
# end
#
# g = zeros(1200)
# inds = SVector{12}([i for i=1:12])
# sv = SVector{12}(zeros(12))
# N = 7
# din = solver_directgames.dyn_constraints.constraints[1].vals
# b = Bundle2(g, din, inds, sv, N)
#
# function dumb2(b::Bundle2{T}) where T
# 	dumb2(b.g, b.din, b.inds, b.sv, b.N)
# 	return nothing
# end
#
#
# function dumb2(g::Vector{T}, din::Vector{SVector{12,T}}, inds::SVector{12,Int}, sv::SVector{12,T}, N::Int) where T
# 	for i=1:N
# 		inds2 = inds .+1
# 		g[inds2.+i] = din[i]
# 	end
# 	return nothing
# end
#
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(g, din, inds, sv, N)
# @allocated  dumb2(b)
# @allocated  dumb2(b)
# @allocated  dumb2(b)
# @allocated  dumb2(b)
# @allocated  dumb2(b)
#
# @btime dumb2(b)
#
#
# @btime $solver_directgames.g_[$solver_directgames.Z[$k]._x] *=2
#
# dyn = solver_directgames.dyn_constraints.constraints[1].vals
# dyn[1]
#
#
#
#
#
#
#
#
#
# pus = [SVector{2}([1,2]),SVector{2}([3,4]),SVector{2}([5,6])]
#
# k = 1
# i = 1
# ix = solver_directgames.Z[k]._x
# iu = solver_directgames.Z[k]._u
# @btime $(solver_directgames.∇F)[$k][$ix,$iu][:,$pus[$i]] # Bki [n,mi]
# fdu = solver_directgames.∇F[k][ix,iu][:,pu[i]] # Bki [n,mi]
#
# @btime $solver_directgames.g_[$solver_directgames.uinds_p[$i][$k]] =
# 	$solver_directgames.C[$i].u[$k] + fdu'*$solver_directgames.ν[$i][$k]
#
#
#
#
#
#
#
# @allocated update_g_lite!(solver_directgames)
# @allocated update_g_lite!(solver_directgames)
# @time update_g_lite!(solver_directgames)
# # @time update_H_!(solver_directgames)
#
# @time timing_solve(solver_directgames)




























#
#
# function update_H_!(solver::DirectGamesSolver, In::SMatrix{ln,ln,T,lnn}) where{ln,lnn,T}
# 	# println("************************HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
#
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
# 	end
#
# 	for con in solver.constraints.constraints
# 		update_H_con!(solver,con,n,m,pu,p,N)
# 	end
#
# 	# Upper-right block
# 	for i = 1:p
# 		spui = solver.spu[i]
# 		for k = 2:N-1
# 			ix = Z[k]._x #ok
# 			fdx = ∇F[k][ix,ix] #ok
# 			H_[solver.xinds_p[i][k], solver.νinds[i][k]] .= fdx' #ok
# 			H_[solver.xinds_p[i][k], solver.νinds[i][k-1]] .= - In #ok
# 			# H_[solver.xinds_p[i][k], solver.νinds[i][k-1]] .= - Matrix{T}(I,n,n) #ok
# 			# H_[solver.xinds_p[i][k], solver.νinds[i][k-1]] .= -Diagonal(In) #ok
# 			# H_[solver.xinds_p[i][k], solver.νinds[i][k-1]] .= -Diagonal(@SVector ones(n)) #ok
# 		end
# 		# H_[solver.xinds_p[i][N], solver.νinds[i][N-1]] .= -Diagonal(@SVector ones(n)) #ok
# 		H_[solver.xinds_p[i][N], solver.νinds[i][N-1]] .= -In #ok
#
# 		for k = 1:N-1
# 			ix = Z[k]._x #ok
# 			iu = Z[k]._u #ok
# 			fdu = ∇F[k][ix,iu][:,spui] #ok
# 			H_[solver.uinds_p[i][k], solver.νinds[i][k]] .= fdu' #ok
# 		end
# 	end
#
# 	# Lower-left block
# 	for k = 2:N-1
# 		ix = Z[k]._x #ok
# 		fdx = ∇F[k][ix,ix] #ok
# 		H_[solver.νinds_p[k-1], solver.xinds[k]] .= -In #ok
# 		H_[solver.νinds_p[k], solver.xinds[k]] .= fdx #ok
# 	end
# 	H_[solver.νinds_p[N-1], solver.xinds[N]] .= -In #ok
#
# 	for k = 1:N-1
# 		ix = Z[k]._x #ok
# 		iu = Z[k]._u #ok
# 		fdu = ∇F[k][ix,iu] #ok
# 		H_[solver.νinds_p[k], solver.uinds[k]] .= fdu #okkk
# 	end
# 	return nothing
# end
