export
    timing_solve

function timing_solve(solver::DirectGamesSolver{T}) where T
    reset!(solver, reset_type=:full)
    solve!(solver)
end

function timing_solve(solver::iLQGamesSolver{T}) where T
    reset!(solver, reset_type=:full)
    solve!(solver)
end

function timing_solve(solver::PenaltyiLQGamesSolver{T}) where T
    reset!(solver, reset_type=:full)
    solve!(solver)
end
