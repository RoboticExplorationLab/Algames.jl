################################################################################
# NewtonCore
################################################################################

mutable struct NewtonCore{Vr,Vrt,SMj,SVhi,SVvi,SAr,SAj,Sd}
	res::Vr               # residual vector
	res_tmp::Vrt          # holds a temporary copy of the residual vector
	jac::SMj              # residual sparse jacobian
	probsize::ProblemSize # size of the problem
	horiz_inds::SVhi      # indices for each variable in an horizontal block
	verti_inds::SVvi      # indices for each variable in an vertical block
	res_sub::SAr          # Residual views dictionary
	jac_sub::SAj          # Jacobian views dictionary
	dyn::Sd               # Dictionary of the Jacobian dynamics indices
end

function NewtonCore(probsize::ProblemSize)
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	S = probsize.S
	res = zeros(S)
	res_tmp = deepcopy(res)
	jac = spzeros(S,S)
	verti_inds = vertical_indices(probsize)
	horiz_inds = horizontal_indices(probsize)
	res_sub = residual_views(res, probsize, verti_inds)
	jac_sub = jacobian_views(jac, probsize, verti_inds, horiz_inds)
	dyn = dynamics_indices(probsize)
	TYPE = typeof.((res, res_tmp, jac, horiz_inds, verti_inds, res_sub, jac_sub, dyn))
	return NewtonCore{TYPE...}(res, res_tmp, jac, probsize,
		horiz_inds, verti_inds, res_sub, jac_sub, dyn)
end

################################################################################
# Indices
################################################################################

function vertical_indices(probsize::ProblemSize)
	N = probsize.N
	n = probsize.n
	p = probsize.p
	mi = probsize.mi
	verti_inds = Dict()
	off = 0
	for i = 1:p
		for k = 1:N-1
			stamp = stampify(:opt, i, :x, 1, k+1)
			verti_inds[stamp] = SVector{n,Int}(off .+ (1:n))
			off += n
			stamp = stampify(:opt, i, :u, i, k)
			verti_inds[stamp] = SVector{mi[i],Int}(off .+ (1:mi[i]))
			off += mi[i]
		end
	end
	for k = 1:N-1
		stamp = stampify(:dyn, 1, :x, 1, k)
		verti_inds[stamp] = SVector{n,Int}(off .+ (1:n))
		off += n
	end
	return verti_inds
end

function horizontal_indices(probsize::ProblemSize)
	N = probsize.N
	n = probsize.n
	p = probsize.p
	mi = probsize.mi

	horiz_inds = Dict()
	off = 0
	for k = 1:N-1
		stamp = stampify(:x, 1, k+1)
		horiz_inds[stamp] = SVector{n,Int}(off .+ (1:n))
		off += n
		for i = 1:p
			stamp = stampify(:u, i, k)
			horiz_inds[stamp] = SVector{mi[i],Int}(off .+ (1:mi[i]))
			off += mi[i]
		end
		for i = 1:p
			stamp = stampify(:λ, i, k)
			horiz_inds[stamp] = SVector{n,Int}(off .+ (1:n))
			off += n
		end
	end
	return horiz_inds
end

function idx(core::NewtonCore, stamp::Stamp)
	vstamp = stampify(stamp.prob, stamp.i0, stamp.n1, stamp.i1, stamp.v1)
	hstamp = stampify(stamp.n2, stamp.i2, stamp.v2)
	verti = vertical_idx(core, vstamp)
	horiz = horizontal_idx(core, hstamp)
	return (verti, horiz)
end

function idx(verti_inds::Dict, horiz_inds::Dict, stamp::Stamp)
	vstamp = stampify(stamp.prob, stamp.i0, stamp.n1, stamp.i1, stamp.v1)
	hstamp = stampify(stamp.n2, stamp.i2, stamp.v2)
	verti = verti_inds[vstamp]
	horiz = horiz_inds[hstamp]
	return (verti, horiz)
end

function horizontal_idx(core::NewtonCore, stamp::HStamp)
	return core.horiz_inds[stamp]
end

function horizontal_idx(horiz_inds::Dict, stamp::HStamp)
	return horiz_inds[stamp]
end

function vertical_idx(core::NewtonCore, stamp::VStamp)
	return core.verti_inds[stamp]
end

function vertical_idx(verti_inds::Dict, stamp::VStamp)
	return verti_inds[stamp]
end


################################################################################
# SubArrays
################################################################################

function residual_views(res::SV, probsize::ProblemSize, verti_inds::Dict) where {SV}
	N = probsize.N
	p = probsize.p

	# Residual Views
	res_sub = Dict{VStamp,SubArray}()
	for prob ∈ (:dyn, :opt)
		for i0 = 1:p
			for n1 in (:x,:u)
				for i1 = 1:p
					for v1 = 1:N
						stamp = stampify(prob, i0, n1, i1, v1)
						if valid(stamp,N,p)
							res_sub[stamp] = view(res, vertical_idx(verti_inds, stamp))
						end
					end
				end
			end
		end
	end
	return res_sub
end

function jacobian_views(jac::SM, probsize::ProblemSize, verti_inds::Dict, horiz_inds::Dict) where {SM}
	N = probsize.N
	p = probsize.p

	# Jacobian Views
	jac_sub = Dict{Stamp,SubArray}()
	for prob ∈ (:dyn, :opt)
		for i0 = 1:p
			for n1 in (:x,:u)
				for i1 = 1:p
					for v1 = 1:N
						for n2 in (:x, :u, :λ)
							for i2 = 1:p
								for v2 = 1:N
									stamp = stampify(prob, i0, n1, i1, v1, n2, i2, v2)
									if valid(stamp,N,p)
										jac_sub[stamp] = view(jac, idx(verti_inds, horiz_inds, stamp)...)
									end
								end
							end
						end
					end
				end
			end
		end
	end
	return jac_sub
end


################################################################################
# Dynamics Indices
################################################################################

function dynamics_indices(probsize::ProblemSize)
	n = probsize.n
	p = probsize.p
	mi = probsize.mi
	pu = probsize.pu
	dyn = Dict()
	dyn[:x] = Dict()
	dyn[:u] = Dict()
	dyn[:x][1] = SVector{n,Int}(1:n)
	for i = 1:p
		dyn[:u][i] = SVector{mi[i],Int}(n .+ pu[i])
	end
	return dyn
end


################################################################################
# Player specific mask
################################################################################

function vertical_mask(core::NewtonCore, i::Int)
	vertical_mask(core.probsize, core.verti_inds, i)
end

function vertical_mask(probsize::ProblemSize, verti_inds::Dict, i::Int)
	N = probsize.N
	p = probsize.p
	# pi = probsize.pz[i]
	pi = Vector(1:probsize.n)

	msk = Vector{Int}()
	stamp = VStamp()

	# Select the optimality indices
	prob = :opt
	i0 = i
	for n1 in (:x,:u)
		for i1 = 1:p
			for v1 = 1:N
				stampify!(stamp, prob, i0, n1, i1, v1)
				if valid(stamp,N,p)
					ind = Vector(vertical_idx(verti_inds, stamp))
					(n1 == :x) && (ind = ind[pi])
					push!(msk, ind...)
				end
			end
		end
	end

	# Select the dynamics indices
	prob = :dyn
	i0 = 1
	n1 = :x
	i1 = 1
	for v1 = 1:N
		stampify!(stamp, prob, i0, n1, i1, v1)
		if valid(stamp,N,p)
			ind = Vector(vertical_idx(verti_inds, stamp)[pi])
			push!(msk, ind...)
		end
	end
	return msk
end


function horizontal_mask(core::NewtonCore, i::Int)
	horizontal_mask(core.probsize, core.horiz_inds, i)
end

function horizontal_mask(probsize::ProblemSize, horiz_inds::Dict, i::Int)
	N = probsize.N
	p = probsize.p
	# pi = probsize.pz[i]
	pi = Vector(1:probsize.n)

	msk = Vector{Int}()
	stamp = HStamp()

	# Select the state indices
	n2 = :x
	i2 = 1
	for v2 = 1:N
		stampify!(stamp, n2, i2, v2)
		if valid(stamp,N,p)
			ind = Vector(horizontal_idx(horiz_inds, stamp)[pi])
			push!(msk, ind...)
		end
	end

	# Select the control indices
	n2 = :u
	i2 = i
	for v2 = 1:N
		stampify!(stamp, n2, i2, v2)
		if valid(stamp,N,p)
			ind = Vector(horizontal_idx(horiz_inds, stamp))
			push!(msk, ind...)
		end
	end

	# Select the multipliers indices
	n2 = :λ
	i2 = i
	for v2 = 1:N
		stampify!(stamp, n2, i2, v2)
		if valid(stamp,N,p)
			ind = Vector(horizontal_idx(horiz_inds, stamp)[pi])
			push!(msk, ind...)
		end
	end
	return msk
end
