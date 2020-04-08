export
    set_penalty!


function set_penalty!(conSet::ConstraintSet{T}, μ::T) where T
	nc = length(conSet.constraints)
	pen = μ*ones(nc)
	set_penalty!(conSet, pen)
    return nothing
end

function set_penalty!(conSet::ConstraintSet{T}, pen::Vector{T}) where T
    for (j,con) in enumerate(conSet.constraints)
        ϕ = con.params.ϕ
        for k in eachindex(con.inds)
            new_μ = SVector{length(con.μ[k])}(pen[j]*ones(length(con.μ[k])))
			# con.μ[k] = clamp.(ϕ * new_μ, 0.0, con.params.μ_max) # used for RSS2020
			con.μ[k] = clamp.(new_μ, 0.0, con.params.μ_max) # used for PFpaper
        end
    end
    return nothing
end

function set_penalty!(solver::PenaltyiLQGamesSolver{T}, pen::Vector{T}) where T
    solver.pen .= pen
    set_penalty!(solver.constraints, pen)
    return nothing
end

function set_penalty!(solver::PenaltyiLQGamesSolver{T}, μ::T) where T
	solver.opts.μ_penalty = μ
    set_penalty!(solver.constraints, μ)
    return nothing
end

function set_penalty!(solver::DirectGamesSolver{T}, pen::Vector{T}) where T
	solver.opts.μ_penalty = pen[1]
    set_penalty!(solver.penalty_constraints, pen)
    return nothing
end

function set_penalty!(solver::DirectGamesSolver{T}, μ::T) where T
	solver.opts.μ_penalty = μ
    set_penalty!(solver.penalty_constraints, μ)
    return nothing
end

function set_penalty!(solver::DirectGamesSolver{T}) where T
    set_penalty!(solver.penalty_constraints, solver.opts.μ_penalty)
    return nothing
end
