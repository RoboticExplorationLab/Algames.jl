# Sample Documentation
```@meta
CurrentModule = Algames
```

This page shows how to pull in the doc strings in your code.

The first thing to do is insert
````markdown
```@meta
CurrentModule = Algames  # your package name here
```
````

which sets the module to your package so you don't have to prepend the methods with your
package.

## Pulling in some docstrings
First, let's pull in the docstrings for `vec_add!` and `vec_sub!`, which we can do using
````markdown
```@docs
vec_add!
vec_sub!
```
````

This inserts the following into our markdown file:

```@docs
vec_add!
vec_sub!
```

Notice how in the `vec_add!` docstring we included both signatures in a single docstring,
but `vec_sub!` had two separate docstrings. We can select only one of the docstrings by
filtering with the input signature:
````markdown
```@docs
vec_sub!(::VecPair)
norm(::VecPair)
```
````

which inserts only one docstring,
```@docs
vec_sub!(::VecPair)
```

## Linking Docstrings
We can link to the docstring for [`vec_add!`](@ref) using the `[`vec_add!`](@ref)` syntax.
Note the tick marks around the method, inside the square brackets. We can also do this
inside the docstring themselves, like we do in the docstring for [`VecPair`](@ref):

```@docs
VecPair
```

For illustration, we also show in this docstring how to include ``\LaTeX`` math inside the
docstring. For reference, we've copied the raw docstring below:

```julia
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
```
