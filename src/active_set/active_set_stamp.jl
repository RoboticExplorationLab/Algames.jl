################################################################################
# Stamps
################################################################################

################################################################################
# Type Definition
################################################################################

abstract type AbstractStamp
end

mutable struct CStamp <: AbstractStamp
	dim::Symbol # dimension (:v or :h)
	con::Symbol # name of the constraint
	i::Int      # player i index
	j::Int      # player j index
	k::Int      # time index k
end

################################################################################
# Base functions
################################################################################

import Base.hash
import Base.==
import Base.isequal

Base.hash(stamp::CStamp, h::UInt)  = hash(stamp.dim, hash(stamp.con, hash(stamp.i, hash(stamp.j, hash(stamp.k, h)))))

function (==)(stamp1::CStamp, stamp2::CStamp)
    out = true
	for name in CStamp.name.names
        out &= getfield(stamp1, name) == getfield(stamp2, name)
    end
    return out
end

Base.isequal(stamp1::CStamp, stamp2::CStamp) = stamp1 == stamp2

################################################################################
# Methods
################################################################################

function CStamp()
	return CStamp(:x, :x, 0, 0, 0)
end

function stampify(dim::Symbol, con::Symbol, i::Int, j::Int, k::Int)
	return CStamp(dim,con,i,j,k)
end

function stampify!(stamp::CStamp, dim::Symbol, con::Symbol, i::Int, j::Int, k::Int)
	stamp.dim = dim
	stamp.con = con
	stamp.i = i
	stamp.j = j
	stamp.k = k
	return nothing
end

################################################################################
# Stamp Validity
################################################################################

function valid(s::CStamp, N::Int, p::Int)
	b = false
	if s.dim == :v # Vertical
		if s.i < s.j
			if s.i ∈ (1:p) && s.j ∈ (1:p) && s.k ∈ (2:N)
				return true
			end
		end
	elseif s.dim == :h # Horizontal
		if s.i ∈ (1:p) && s.j ∈ (1:p) && s.k ∈ (2:N) && s.i != s.j
			return true
		end
	else
		@show "Error invalid dim Symbol."
	end
	return b
end
