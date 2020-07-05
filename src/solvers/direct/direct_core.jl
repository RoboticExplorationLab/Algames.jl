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
    for c in eachindex(conSet.constraints)
        TO.dual_update!(conSet.constraints[c])
    end
end

function penalty_update!(solver::DirectGamesSolver)
	# @show "penalty update"
    for c in eachindex(solver.constraints.constraints)
        TO.penalty_update!(solver.constraints.constraints[c])
    end
	for c in eachindex(solver.dyn_constraints.constraints)
		TO.penalty_update!(solver.dyn_constraints.constraints[c])
	end
end

function inner_step!(solver::DirectGamesSolver)
	# @show "inner step"
	solver.δY .= -(lu(solver.H_)\solver.g_) # better
	# solver.δY .= -solver.H_\solver.g_
	α, optimality_merit = line_search!(solver)
	record_inner_iteration!(solver, α)
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

		###
		###
		TO.evaluate!(solver.penalty_constraints, solver.Z̄)
		TO.jacobian!(solver.penalty_constraints, solver.Z̄)
		TO.update_active_set!(solver.penalty_constraints,solver.Z̄)
		penalty_expansion!(solver, solver.Z̄)
		###
		###

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
	p = get_model(solver).p
	pu = get_model(solver).pu
	N = solver.N
	for i = 1:p
        cost!(solver.obj[i], solver.Z, fill(Array(pu[i]),N))
    end
	_J = TO.get_J.(solver.obj)
	J = sum.(_J)
	return J
end

# function penalty_expansion!(solver::DirectGamesSolver{T}, Z::TO.Traj) where T
#     C = solver.C
#     conSet = solver.penalty_constraints
#     n,m,pu,p = size(solver.model)
#     n,m,N = size(solver)
#     for i = 1:p
# 		for c in eachindex(conSet.constraints)
# 			con = conSet.constraints[c]
#             if typeof(con).parameters[2] == State
#                 for (j,k) in enumerate(intersect(con.inds,1:N))
#                     Iμ = TO.penalty_matrix(con,j)
#                     C[i].xx[k] += con.∇c[j]'*Iμ*con.∇c[j]
#                     C[i].x[k] += con.∇c[j]'*Iμ*con.vals[j]
#                 end
#             elseif typeof(con).parameters[2] == Control
#                 for (j,k) in enumerate(intersect(con.inds,1:N-1))
#                     Iμ = TO.penalty_matrix(con,j)
#                     C[i].uu[k] += con.∇c[j]'*Iμ*con.∇c[j]
#                     C[i].u[k] += con.∇c[j]'*Iμ*con.vals[j]
#                 end
# 			else
# 				@warn(typeof(con), " this type of constraint is not supported.")
#             end
#         end
#     end
#     return nothing
# end


function penalty_expansion!(solver::DirectGamesSolver{T}, Z::TO.Traj) where T
    C = solver.C
    conSet = solver.penalty_constraints
    n,m,pu,p = size(solver.model)
    n,m,N = size(solver)
	for c in eachindex(conSet.constraints)
		con = conSet.constraints[c]
		update_penalty_expansion!(solver,con,n,m,p,N)
    end
    return nothing
end

function update_penalty_expansion!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,State},
	n::Int, m::Int, p::Int, N::Int) where {T}
	C = solver.C
	for i=1:p
	    for (j,k) in enumerate(intersect(con.inds,1:N))
			Iμ = TO.penalty_matrix(con,j)
			C[i].xx[k] += con.∇c[j]'*Iμ*con.∇c[j]
			C[i].x[k] += con.∇c[j]'*Iμ*con.vals[j]
		end
	end
	return nothing
end

function update_penalty_expansion!(solver::DirectGamesSolver{T}, con::ConstraintVals{T,Control},
	n::Int, m::Int, p::Int, N::Int) where {T}
	C = solver.C
	for i=1:p
	    for (j,k) in enumerate(intersect(con.inds,1:N-1))
	        Iμ = TO.penalty_matrix(con,j)
	        C[i].uu[k] += con.∇c[j]'*Iμ*con.∇c[j]
	        C[i].u[k] += con.∇c[j]'*Iμ*con.vals[j]
	    end
	end
	return nothing
end
