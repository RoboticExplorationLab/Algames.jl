################################################################################
# Stamps
################################################################################

################################################################################
# Type Definition
################################################################################

abstract type AbstractStamp
end

mutable struct Stamp <: AbstractStamp
	prob::Symbol # name of the problem
	i0::Int      # index of the problem
	n1::Symbol   # name of the variable
	i1::Int      # index of the variable
	v1::Int      # value of the variable
	n2::Symbol   # name of the variable
	i2::Int      # index of the variable
	v2::Int      # value of the variable
end

mutable struct VStamp <: AbstractStamp
	prob::Symbol # name of the problem
	i0::Int      # index of the problem
	n1::Symbol   # name of the variable
	i1::Int      # index of the variable
	v1::Int      # value of the variable
end

mutable struct HStamp <: AbstractStamp
	n2::Symbol   # name of the variable
	i2::Int      # index of the variable
	v2::Int      # value of the variable
end

################################################################################
# Base functions
################################################################################

import Base.hash
import Base.==
import Base.isequal

Base.hash(stamp::Stamp, h::UInt)  = hash(stamp.prob, hash(stamp.i0, hash(stamp.n1, hash(stamp.i1, hash(stamp.v1, hash(stamp.n2, hash(stamp.i2, hash(stamp.v2, h))))))))
Base.hash(stamp::VStamp, h::UInt) = hash(stamp.prob, hash(stamp.i0, hash(stamp.n1, hash(stamp.i1, hash(stamp.v1, h)))))
Base.hash(stamp::HStamp, h::UInt) = hash(stamp.n2, hash(stamp.i2, hash(stamp.v2, h)))

function (==)(stamp1::AbstractStamp, stamp2::AbstractStamp)
    out = true
    T1 = typeof(stamp1)
    T2 = typeof(stamp2)
    T1 == T2 ? nothing : return false
    for name in fieldnames(T1)
        out &= getfield(stamp1, name) == getfield(stamp2, name)
    end
    return out
end

function (==)(stamp1::Stamp, stamp2::Stamp)
    out = true
	for name in Stamp.name.names
        out &= getfield(stamp1, name) == getfield(stamp2, name)
    end
    return out
end

function (==)(stamp1::VStamp, stamp2::VStamp)
    out = true
	for name in VStamp.name.names
        out &= getfield(stamp1, name) == getfield(stamp2, name)
    end
    return out
end

function (==)(stamp1::HStamp, stamp2::HStamp)
    out = true
	for name in HStamp.name.names
        out &= getfield(stamp1, name) == getfield(stamp2, name)
    end
    return out
end

Base.isequal(stamp1::Stamp, stamp2::Stamp) = stamp1 == stamp2
Base.isequal(stamp1::VStamp, stamp2::VStamp) = stamp1 == stamp2
Base.isequal(stamp1::HStamp, stamp2::HStamp) = stamp1 == stamp2

################################################################################
# Methods
################################################################################

function Stamp()
	return Stamp(:x, 0, :x, 0, 0, :x, 0, 0)
end

function VStamp()
	return VStamp(:x, 0, :x, 0, 0)
end

function HStamp()
	return HStamp(:x, 0, 0)
end

function stampify(name::Symbol, step::Int)
	if name == :x_1
		name = :x
		step -= 1
	elseif name == :x1
		name = :x
		step += 1
	end
	return name, step
end

function stampify(prob::Symbol, ind_p::Int, name_i::Symbol, ind_i::Int, step_i::Int, name_j::Symbol, ind_j::Int, step_j::Int)
	name_i, step_i = stampify(name_i, step_i)
	name_j, step_j = stampify(name_j, step_j)
	return Stamp(prob, ind_p, name_i, ind_i, step_i, name_j, ind_j, step_j)
end

function stampify!(stamp::Stamp, prob::Symbol, ind_p::Int, name_i::Symbol, ind_i::Int, step_i::Int, name_j::Symbol, ind_j::Int, step_j::Int)
	name_i, step_i = stampify(name_i, step_i)
	name_j, step_j = stampify(name_j, step_j)
	stamp.prob = prob
	stamp.i0 = ind_p
	stamp.n1 = name_i
	stamp.i1 = ind_i
	stamp.v1 = step_i
	stamp.n2 = name_j
	stamp.i2 = ind_j
	stamp.v2 = step_j
	return nothing
end

function stampify(prob::Symbol, ind_p::Int, name_i::Symbol, ind_i::Int, step_i::Int)
	name_i, step_i = stampify(name_i, step_i)
	return VStamp(prob, ind_p, name_i, ind_i, step_i)
end

function stampify!(stamp::VStamp, prob::Symbol, ind_p::Int, name_i::Symbol, ind_i::Int, step_i::Int)
	name_i, step_i = stampify(name_i, step_i)
	stamp.prob = prob
	stamp.i0 = ind_p
	stamp.n1 = name_i
	stamp.i1 = ind_i
	stamp.v1 = step_i
	return nothing
end

function stampify(name_j::Symbol, ind_j::Int, step_j::Int)
	name_j, step_j = stampify(name_j, step_j)
	return HStamp(name_j, ind_j, step_j)
end

function stampify!(stamp::HStamp, name_j::Symbol, ind_j::Int, step_j::Int)
	name_j, step_j = stampify(name_j, step_j)
	stamp.n2 = name_j
	stamp.i2 = ind_j
	stamp.v2 = step_j
	return nothing
end

################################################################################
# Stamp Validity
################################################################################

function valid(s::Stamp, N::Int, p::Int)
	return valid(s.prob, s.i0, s.n1, s.i1, s.v1, s.n2, s.i2, s.v2, N, p)
end

function valid(prob::Symbol, i0::Int, n1::Symbol, i1::Int, v1::Int, n2::Symbol, i2::Int, v2::Int, N::Int, p::Int)
	b1 = valid(prob, i0, n1, i1, v1, N, p) # stamp 1 is valid
	b2 = false # stamp 2 is valid

	if prob == :opt && i0 ∈ (1:p)
		if n2 == :u && i2 ∈ (1:p) && v2 ∈ (1:N-1) # uj1, ...N-1
			b2 = true
		elseif n2 == :λ && i2 == i0 && v2 ∈ (1:N-1) # λi1 ...N-1
			b2 = true
		elseif n2 == :x && i2 == 1 && v2 ∈ (2:N) # x2...xN
			b2 = true
		end
	end
	if prob == :dyn && i0 == 1
		if n2 == :u && i2 ∈ (1:p) && v2 ∈ (1:N-1) # uj1, ...N-1
			b2 = true
		elseif n2 == :x && i2 == 1 && v2 ∈ (2:N) # x2...xN
			b2 = true
		end
	end
	return b1 && b2
end

function valid(s::VStamp, N::Int, p::Int)
	return valid(s.prob, s.i0, s.n1, s.i1, s.v1, N, p)
end

function valid(prob::Symbol, i0::Int, n1::Symbol, i1::Int, v1::Int, N::Int, p::Int)
	b1 = false # stamp 1 is valid
	if prob == :opt && i0 ∈ (1:p)
		if n1 == :u && i1 == i0 && v1 ∈ (1:N-1) # u1...uN-1
			b1 = true
		elseif n1 == :x && i1 == 1 && v1 ∈ (2:N) # x2...xN
			b1 = true
		end
	end
	if prob == :dyn && i0 == 1
		if n1 == :x && i1 == 1 && v1 ∈ (1:N-1) # x1...xN-1
			b1 = true
		end
	end
	return b1
end

function valid(s::HStamp, N::Int, p::Int)
	return valid(s.n2, s.i2, s.v2, N, p)
end

function valid(n2::Symbol, i2::Int, v2::Int, N::Int, p::Int)
	b2 = false # stamp 2 is valid
	if n2 == :u && i2 ∈ (1:p) && v2 ∈ (1:N-1) # uj1, ...N-1
		b2 = true
	elseif n2 == :λ && i2 ∈ (1:p) && v2 ∈ (1:N-1) # λi1 ...N-1
		b2 = true
	elseif n2 == :x && i2 == 1 && v2 ∈ (2:N) # x2...xN
		b2 = true
	end
	return b2
end
