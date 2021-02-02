################################################################################
# Active Set
################################################################################

function active(game_con::GameConstraintValues, stamp::CStamp)
	probsize = game_con.probsize
    N = probsize.N
    p = probsize.p
    px = probsize.px

    a = BitArray([0])
    if stamp.con == :col
        if valid(stamp, N, p)
            i = stamp.i
            j = stamp.j
            j_ = findfirst(x -> all(x.con.x2 .== px[j]), game_con.state_conval[i])
            conval = game_con.state_conval[i][j_]
            l = findfirst(x -> x==stamp.k, conval.inds)
            a[1] = conval.active[l][1]
			# if conval.vals[l][1] >= 0. || conval.λ[l][1] != 0.
                # a[1] = 1
            # end
        end
    end
    return a
end

function active_vertical_mask!(ascore::ActiveSetCore, game_con::GameConstraintValues)
    probsize = game_con.probsize
    N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	S = probsize.S

    ascore.vmask = Vector{Int}(1:S)
	stamp = CStamp()
	for k = 2:N
		for i = 1:p
			for j = i+1:p
				stampify!(stamp, :v, :col, i, j, k)
				if active(game_con, stamp)[1] == 1
					ascore.vmask = vcat(ascore.vmask, vertical_idx(ascore, stamp))
				end
			end
		end
	end
    return nothing
end

function active_horizontal_mask!(ascore::ActiveSetCore, game_con::GameConstraintValues)
    probsize = game_con.probsize
    N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	S = probsize.S

    ascore.hmask = Vector{Int}(1:S)
	stamp = CStamp()
	for k = 2:N
		for i = 1:p
			for j ∈ setdiff(1:p,i)
				stampify!(stamp, :h, :col, i, j, k)
				if active(game_con, stamp)[1] == 1
					ascore.hmask = vcat(ascore.hmask, horizontal_idx(ascore, stamp))
				end
			end
		end
	end
    return nothing
end


################################################################################
# Residual
################################################################################

function get_collision_conval(game_con::GameConstraintValues, i::Int, j::Int)
	px = game_con.probsize.px
	j_ = findfirst(
		x -> typeof(x.con) <: CollisionConstraint && all(x.con.x2 .== px[j]),
		game_con.state_conval[i]
	)
	if j_ == nothing
		valid = false
		return nothing, valid
	else
		valid = true
		return game_con.state_conval[i][j_], valid
	end
end

function residual!(ascore::ActiveSetCore, prob::GameProblem{KN,n,m,T,SVd,SVx},
	pdtraj::PrimalDualTraj{KN,n,m,T,SVd}) where {KN,n,m,T,SVd,SVx}

	probsize = ascore.probsize
	N = probsize.N
	p = probsize.p
	px = probsize.px
	S = probsize.S

	# Initialization
	ascore.res .= 0.0

	residual!(prob, pdtraj)
	ascore.res[1:S] .= prob.core.res

	evaluate!(prob.game_con, pdtraj.pr)
	stamp = CStamp()
	for i = 1:p
		for j = i+1:p
			conval, v = get_collision_conval(prob.game_con, i, j)
			if v
				for (l,k) in enumerate(conval.inds)
					stampify!(stamp, :v, :col, i, j, k)
					valid(stamp,N,p) ? add2sub(ascore.res_sub[stamp], conval.vals[l]) : nothing
				end
			end
		end
	end
    return nothing
end


################################################################################
# Residual Jacobian
################################################################################

function residual_jacobian!(ascore::ActiveSetCore, prob::GameProblem{KN,n,m,T,SVd,SVx},
	pdtraj::PrimalDualTraj{KN,n,m,T,SVd}) where {KN,n,m,T,SVd,SVx}

	probsize = ascore.probsize
	N = probsize.N
	p = probsize.p
	S = probsize.S

	# Reset!
	sparse_zero!(ascore.jac)

	residual_jacobian!(prob, pdtraj)
	ascore.jac[1:S,1:S] .= prob.core.jac

    # Allocations
	vs = VStamp()
	hs = HStamp()
    cs = CStamp()

	jacobian!(prob.game_con, pdtraj.pr)
	for i = 1:p
		for j ∈ setdiff(1:p,i)
			conval, v = get_collision_conval(prob.game_con, i, j)
			if v
				for (l,k) in enumerate(conval.inds)
					stampify!(vs, :opt, i, :x, 1, k)
					stampify!(cs, :h, :col, i, j, k)
					if valid(cs,N,p) && valid(vs,N,p)
						add2sub(ascore.jac_sub[(vs,cs)], conval.jac[l]')
					end
					if valid(hs,N,p) && valid(cs,N,p)
						stampify!(hs, :x, 1, k)
						stampify!(cs, :v, :col, i, j, k)
						add2sub(ascore.jac_sub[(hs,cs)], conval.jac[l])
					end
				end
			end
		end
	end
    return nothing
end


function update_nullspace!(ascore::ActiveSetCore, prob::GameProblem{KN,n,m,T,SVd,SVx},
	pdtraj::PrimalDualTraj{KN,n,m,T,SVd}; atol::T=1e-20) where {KN,n,m,T,SVd,SVx}

	update_active_set!(prob.game_con, pdtraj.pr)
	active_vertical_mask!(ascore, prob.game_con)
	active_horizontal_mask!(ascore, prob.game_con)
	residual_jacobian!(ascore, prob, pdtraj)
	djac = Matrix(ascore.jac[ascore.vmask, ascore.hmask])
	reset!(ascore.null)
	add_matrix!(ascore.null, nullspace(djac, atol=1e-20), ascore.hmask)
	return nothing
end
