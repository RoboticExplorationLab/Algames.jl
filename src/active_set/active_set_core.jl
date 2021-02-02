################################################################################
# NullSpace
################################################################################

mutable struct NullSpace{T,SVv,SVΔ,SVλ}
	probsize::ProblemSize
	mat::Matrix{T}
	vec::Vector{SVv}
	Δtraj::Vector{SVΔ}
	Δλ::Vector{SVλ}
end

function NullSpace(probsize)
	mat = zeros(0,0)
	vec = Vector([zeros(0)])
	Δtraj = Vector([zeros(0)])
	Δλ = Vector([zeros(0)])
	return NullSpace(probsize, mat, vec, Δtraj, Δλ)
end

function reset!(null::NullSpace{T,SVv,SVΔ,SVλ}) where {T,SVv,SVΔ,SVλ}
	null.mat = zeros(T,0,0)
	null.vec = Vector{SVv}()
	null.Δtraj = Vector{SVΔ}()
	null.Δλ = Vector{SVλ}()
	return nothing
end

function add_matrix!(null::NullSpace, mat::AbstractMatrix, hmask::Vector{Int})
	probsize = null.probsize
	N = probsize.N
	p = probsize.p
	S = probsize.S
	Sh = S + p*(p-1)*(N-1)
	m1, m2 = size(mat)
	@assert (S <= m1 <= Sh)

	null.mat = mat
	for l = 1:m2
		vec = zeros(Sh)
		vec[hmask] = mat[:,l]
		vec ./= mean(abs.(vec))
		Δtraj = vec[1:S]
		Δλ = vec[S+1:end]
		push!(null.vec, vec)
		push!(null.Δtraj, Δtraj)
		push!(null.Δλ, Δλ)
	end
	return nothing
end

################################################################################
# ActiveSetCore
################################################################################

mutable struct ActiveSetCore{Vr,Vrt,SMj,SVhi,SVvi,SAr,SAj,Sd}
	res::Vr               # residual vector
	res_tmp::Vrt          # holds a temporary copy of the residual vector
	jac::SMj              # residual sparse jacobian
	probsize::ProblemSize # size of the problem
	horiz_inds::SVhi      # indices for each variable in an horizontal block
	verti_inds::SVvi      # indices for each variable in an vertical block
	res_sub::SAr          # Residual views dictionary
	jac_sub::SAj          # Jacobian views dictionary
	dyn::Sd               # Dictionary of the Jacobian dynamics indices
	vmask::Vector{Int}    # Vertical mask of the active inequality constraints
	hmask::Vector{Int}    # Horizontal mask of the active inequality constraints
	null::NullSpace       # Structure holding nullspace of the active set Jacobian
end

function ActiveSetCore(probsize::ProblemSize)
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	Sv = probsize.S + Int(p*(p-1)*(N-1)/2)
	Sh = probsize.S + p*(p-1)*(N-1)
	res = zeros(Sv)
	res_tmp = deepcopy(res)
	jac = spzeros(Sv,Sh)
	verti_inds = complete_vertical_indices(probsize)
	horiz_inds = complete_horizontal_indices(probsize)
	res_sub = complete_residual_views(res, probsize, verti_inds)
	jac_sub = complete_jacobian_views(jac, probsize, verti_inds, horiz_inds)
	dyn = dynamics_indices(probsize)
	vmask = Vector{Int}(1:Sv)
	hmask = Vector{Int}(1:Sh)
	null = NullSpace(probsize)
	TYPE = typeof.((res, res_tmp, jac, horiz_inds, verti_inds, res_sub, jac_sub, dyn))
	return ActiveSetCore{TYPE...}(res, res_tmp, jac, probsize,
		horiz_inds, verti_inds, res_sub, jac_sub, dyn, vmask, hmask, null)
end

################################################################################
# Indices
################################################################################

function complete_vertical_indices(probsize::ProblemSize)
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
	for k = 2:N
		for i = 1:p
			for j = i+1:p
				stamp = stampify(:v, :col, i, j, k)
				verti_inds[stamp] = SVector{1,Int}(off .+ (1:1))
				off += 1
			end
		end
	end
	return verti_inds
end

function complete_horizontal_indices(probsize::ProblemSize)
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
	for k = 2:N
		for i = 1:p
			for j ∈ setdiff(1:p,i)
				stamp = stampify(:h, :col, i, j, k)
				horiz_inds[stamp] = SVector{1,Int}(off .+ (1:1))
				off += 1
			end
		end
	end
	return horiz_inds
end

function horizontal_idx(core::ActiveSetCore, stamp::CStamp)
	return core.horiz_inds[stamp]
end

function horizontal_idx(horiz_inds::Dict, stamp::CStamp)
	return horiz_inds[stamp]
end

function vertical_idx(core::ActiveSetCore, stamp::CStamp)
	return core.verti_inds[stamp]
end

function vertical_idx(verti_inds::Dict, stamp::CStamp)
	return verti_inds[stamp]
end


################################################################################
# SubArrays
################################################################################

function complete_residual_views(res::SV, probsize::ProblemSize, verti_inds::Dict) where {SV}
	N = probsize.N
	p = probsize.p

	# Residual Views
	res_sub = Dict{AbstractStamp,SubArray}()
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
	for k = 2:N
		for i = 1:p
			for j = i+1:p
				stamp = stampify(:v,:col, i, j, k)
				if valid(stamp,N,p)
					res_sub[stamp] = view(res, vertical_idx(verti_inds, stamp))
				end
			end
		end
	end
	return res_sub
end

function complete_jacobian_views(jac::SM, probsize::ProblemSize, verti_inds::Dict, horiz_inds::Dict) where {SM}
	N = probsize.N
	p = probsize.p

	# Jacobian Views
	jac_sub = Dict{Any,SubArray}()

	for prob ∈ (:opt,)
		for i0 = 1:p
			for n1 in (:x,)
				for i1 = 1:1
					for v1 = 1:N
						for j = 1:p
							vstamp = stampify(prob, i0, n1, i1, v1)
							hstamp = stampify(:h, :col, i0, j, v1)
							if valid(vstamp,N,p) && valid(hstamp,N,p)
								jac_sub[(vstamp, hstamp)] = view(jac,
									vertical_idx(verti_inds, vstamp),
									horizontal_idx(horiz_inds, hstamp))
							end
						end
					end
				end
			end
		end
	end
	for n2 in (:x,)
		for i2 = 1:1
			for v2 = 2:N
				for j = i2+1:p
					vstamp = stampify(:v, :col, i2, j, v2)
					hstamp = stampify(n2, i2, v2)
					if valid(vstamp,N,p) && valid(hstamp,N,p)
						jac_sub[(vstamp, hstamp)] = view(jac,
							vertical_idx(verti_inds, vstamp),
							horizontal_idx(horiz_inds, hstamp))
					end
				end
			end
		end
	end
	return jac_sub
end
