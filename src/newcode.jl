"""
    vec_add!(a::AbstractVector, b::AbstractVector)
    ved_add!(v::VecPair)

Adds `b` to `a`, modifying `a` in place. Returns the modified `a`.
Vectors must be the same length.
"""
function vec_add!(a::AbstractVector, b::AbstractVector)
    @assert length(a) == length(b)
    for i in eachindex(a)
        a[i] += b[i]
    end
    return a
end

"""
    vec_sub!(a::AbstractVector, b::AbstractVector)

Substracts `b` from `a`, modifying `a` in place. Returns the modified `a`.
Vectors must be the same length.
"""
function vec_sub!(a::AbstractVector, b::AbstractVector)
    @assert length(a) == length(b)
    for i in eachindex(a)
        a[i] -= b[i]
    end
    return a
end

"""
    VecPair{V}

Holds two vectors of the same length and type.

The vectors can be retrieved using `v.a` and `v.b` or `v[1]` and `v[2]`.
Supports [`vec_add!`](@ref) and [`vec_sub!`](@ref).

Here is some ``\\LaTeX`` for you:
```math
    \\sum_{i=1}^N x_k^T Q_k x_k
```

# Constructors
    VecPair{V}(a,b)
    VecPair(a::V, b::V)
    VecPair(a::StaticVector, b::StaticVector)

"""
struct VecPair{V}
    a::V
    b::V
    function VecPair{V}(a,b) where V <: AbstractVector
        @assert length(a) == length(b)
        new{V}(a,b)
    end
end

VecPair(a::V, b::V) where V <: AbstractVector = VecPair{V}(a,b)
VecPair(a::StaticVector{N,T}, b::StaticVector{N,T}) where {N,T} = VecPair{MVector{N,T}}(a,b)

vec_add!(v::VecPair) = vec_add!(v.a, v.b)

"""
    vec_sub!(v::VecPair)

Substract `v.b` from `v.a`, modifying `v.a` in place.
"""
vec_sub!(v::VecPair) = vec_sub!(v.a, v.b)

"""
    norm(v::VecPair, p=2)

Return the `p`-norm of `[v.a; v.b]`.
"""
LinearAlgebra.norm(v::VecPair, p::Real=2) = norm(SA[norm(v.a, p), norm(v.b, p)], p)
function Base.getindex(v::VecPair, i::Int)
    if i == 1
        v.a
    elseif i == 2
        v.b
    else
        throw(ArgumentError("Index can only be 1 or 2, got $i"))
    end
end
