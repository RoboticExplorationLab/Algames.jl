export
	solve!,
	step!,
	record_iteration!,
	set_tolerances!,
	evaluate_convergence,
	dual_update!,
	penalty_update!


function TO.solve!(solver::ALGamesSolver{T,S}) where {T,S}
	TO.set_verbosity!(solver.opts)
	TO.clear_cache!(solver.opts)
    c_max::T = Inf

	# Extract stuff from solver
	Z = get_trajectory(solver)
	obj = TO.get_objective(solver.solver_uncon)#::Vector{ALObjective{T}}
    conSets = [obj_i.constraints for obj_i in obj]
    solver_uncon = solver.solver_uncon::S

	# Reset solver
    TO.reset!.(conSets)
    solver.stats.iterations = 0
	solver.stats.iterations_total = 0

	n,m,pl,p = size(solver.solver_uncon.model)
	n,m,N = size(solver)
	# Calculate cost
    # cost!(obj, Z)
	for i = 1:p
        cost!(obj[i], Z, fill(pl[i],N))
    end
    J_ = TO.get_J.(obj)
    J = sum.(J_)

    for j = 1:solver.opts.iterations
        TO.set_tolerances!(solver,solver_uncon,j)

        step!(solver)
        J = sum.(J_)
		c_max = 0.0
		for i = 1:p
			c_max = max(c_max, maximum(conSets[i].c_max))
		end

        record_iteration!(solver, J, c_max)

        converged = evaluate_convergence(solver)
        if converged
            break
        end

        reset!(solver_uncon, false)
    end
    return solver
end

function step!(solver::ALGamesSolver)

    # Solve the unconstrained problem
    solve!(solver.solver_uncon)

    # Outer loop update
    TO.dual_update!(solver)
    TO.penalty_update!(solver)
    TO.max_violation!(get_constraints(solver))

	# Reset verbosity level after it's modified
	TO.set_verbosity!(solver.opts)
end

function record_iteration!(solver::ALGamesSolver{T,S}, J::Vector{T}, c_max::T) where {T,S}

    solver.stats.iterations += 1
    i = solver.stats.iterations
    solver.stats.iterations_total += solver.solver_uncon.stats.iterations
    solver.stats.c_max[i] = c_max
	solver.stats.cost[i] = J

	conSet = TO.get_constraints(solver)
	TO.max_penalty!(conSet)
	solver.stats.penalty_max[i] = maximum(conSet.c_max)

	@logmsg OuterLoop :iter value=i
	@logmsg OuterLoop :total value=solver.stats.iterations_total
	@logmsg OuterLoop :cost value=J
    @logmsg OuterLoop :c_max value=c_max
	if solver.opts.verbose
		print_level(OuterLoop)
	end
end

function TO.set_tolerances!(solver::ALGamesSolver{T},
        solver_uncon::TO.AbstractSolver{T},i::Int) where T
    if i != solver.opts.iterations
        solver_uncon.opts.cost_tolerance = solver.opts.cost_tolerance_intermediate
        solver_uncon.opts.gradient_norm_tolerance = solver.opts.gradient_norm_tolerance_intermediate
    else
        solver_uncon.opts.cost_tolerance = solver.opts.cost_tolerance
        solver_uncon.opts.gradient_norm_tolerance = solver.opts.gradient_norm_tolerance
    end

    return nothing
end

function evaluate_convergence(solver::ALGamesSolver)
	i = solver.stats.iterations
    solver.stats.c_max[i] < solver.opts.constraint_tolerance ||
		solver.stats.penalty_max[i] >= solver.opts.penalty_max
end

"General Dual Update"
function dual_update!(solver::ALGamesSolver) where {T,Q,N,M,NM}
    conSet = get_constraints(solver)
    for i in eachindex(conSet.constraints)
        dual_update!(conSet.constraints[i])
    end
end

"Dual Update for Equality Constraints"
function dual_update!(con::ConstraintVals{T,W,C}) where
		{T,W,C<:TO.AbstractConstraint{Equality}}
	λ = con.λ
	c = con.vals
	μ = con.μ
	λ_max = con.params.λ_max
	for i in eachindex(con.inds)
		λ[i] = clamp.(λ[i] + μ[i] .* c[i], -λ_max, λ_max)
	end
end

"Dual Update for Inequality Constraints"
function dual_update!(con::ConstraintVals{T,W,C}) where
		{T,W,C<:TO.AbstractConstraint{Inequality}}
	λ = con.λ
	c = con.vals
	μ = con.μ
	for i in eachindex(con.inds)
		λ[i] = clamp.(λ[i] + μ[i] .* c[i], 0.0, con.params.λ_max)
	end
end

"General Penalty Update"
function penalty_update!(solver::ALGamesSolver)
    conSet = get_constraints(solver)
    for i in eachindex(conSet.constraints)
        penalty_update!(conSet.constraints[i])
    end
end

"Penalty Update for ConstraintVals"
function penalty_update!(con::ConstraintVals{T}) where T
	ϕ = con.params.ϕ
	μ = con.μ
	for i in eachindex(con.inds)
		μ[i] = clamp.(ϕ * μ[i], 0.0, con.params.μ_max)
	end
end











#
# function TO.solve!(prob::StaticALGamesProblem, solver::StaticALGamesSolver{T,S}) where {T,S}
#     c_max::T = Inf
#     conSet = prob.obj.constraints
#     reset!(conSet, solver.opts)
#     solver.stats.iterations = 0
#     solver_uncon = solver.solver_uncon::S
#     cost!(prob.obj, prob.Z)
#     J_ = get_J(prob.obj)
#     J = sum(J_)
#
#
#     for i = 1:solver.opts.iterations
#         set_tolerances!(solver,solver_uncon,i)
#
#         step!(prob, solver)
#         J = sum(J_)
#         c_max = maximum(conSet.c_max)
#
#
#         record_iteration!(prob, solver, J, c_max)
#
#         converged = evaluate_convergence(solver)
#         if converged
#             break
#         end
#
#         reset!(solver_uncon, false)
#     end
#     return solver
# end
#
# function step!(prob::StaticALGamesProblem, solver::StaticALGamesSolver)
#
#     # Solve the unconstrained problem
#     solve!(prob, solver.solver_uncon)
#
#     # Outer loop update
#     dual_update!(prob, solver)
#     penalty_update!(prob, solver)
#     max_violation!(prob.obj.constraints)
#     # copyto!(solver.C_prev,solver.C)
#
# end
#
# function record_iteration!(prob::StaticGameProblem, solver::StaticALGamesSolver{T,S},
#         J::T, c_max::T) where {T,S}
#
#     solver.stats.iterations += 1
#     i = solver.stats.iterations
#     solver.stats.iterations_total += solver.solver_uncon.stats.iterations
#     solver.stats.c_max[i] = c_max
# 	solver.stats.cost[i] = J
#     # solver.stats.
#     # push!(solver.stats[:iterations_inner], unconstrained_solver.stats[:iterations])
#     # push!(solver.stats[:cost],J)
#     # push!(solver.stats[:c_max],c_max)
#     # push!(solver.stats[:penalty_max],max_penalty(solver))
#     # push!(solver.stats_uncon, copy(unconstrained_solver.stats))
#     #
# end
#
# function set_tolerances!(solver::StaticALGamesSolver{T},
#         solver_uncon::TO.AbstractSolver{T},i::Int) where T
#     if i != solver.opts.iterations
#         solver_uncon.opts.cost_tolerance = solver.opts.cost_tolerance_intermediate
#         solver_uncon.opts.gradient_norm_tolerance = solver.opts.gradient_norm_tolerance_intermediate
#     else
#         solver_uncon.opts.cost_tolerance = solver.opts.cost_tolerance
#         solver_uncon.opts.gradient_norm_tolerance = solver.opts.gradient_norm_tolerance
#     end
#
#     return nothing
# end
#
# function evaluate_convergence(solver::StaticALGamesSolver)
#     solver.stats.c_max[solver.stats.iterations] < solver.opts.constraint_tolerance
# end
#
# function dual_update!(con::TO.ConstraintVals{T,W,C},
# 		opts::StaticALGamesSolverOptions{T}) where
# 		{T,W,C<:TO.AbstractStaticConstraint{Equality}}
# 	λ = con.λ
# 	c = con.vals
# 	μ = con.μ
# 	for i in eachindex(con.inds)
# 		λ[i] = clamp.(λ[i] + μ[i] .* c[i], -opts.dual_max, opts.dual_max)
# 	end
# end
#
# function dual_update!(prob::StaticALGamesProblem,
#         solver::StaticALSolver) where {T,Q,N,M,NM}
#     conSet = prob.obj.constraints
#     for i in eachindex(conSet.constraints)
#         dual_update!(conSet.constraints[i], solver.opts)
#     end
# end
#
# function dual_update!(con::TO.ConstraintVals{T,W,C},
# 		opts::StaticALGamesSolverOptions{T}) where
# 		{T,W,C<:TO.AbstractStaticConstraint{Inequality}}
# 	λ = con.λ
# 	c = con.vals
# 	μ = con.μ
# 	for i in eachindex(con.inds)
# 		λ[i] = clamp.(λ[i] + μ[i] .* c[i], 0.0, opts.dual_max)
# 	end
# end
#
#
# function penalty_update!(prob::StaticALGamesProblem, solver::StaticALGamesSolver)
#     conSet = prob.obj.constraints
#     for i in eachindex(conSet.constraints)
#         penalty_update!(conSet.constraints[i], solver.opts)
#     end
# end
#
#
# function penalty_update!(con::ConstraintVals{T}, opts::StaticALGamesSolverOptions{T}) where T
# 	ϕ = opts.penalty_scaling
# 	μ = con.μ
# 	for i in eachindex(con.inds)
# 		μ[i] = clamp.(ϕ * μ[i], 0.0, opts.penalty_max)
# 	end
# end
#
# function max_violation(prob::StaticALGamesProblem)
#     conSet = prob.obj.constraints
#     evaluate!(conSet, prob.Z)
#     max_violation!(conSet)
#     return maximum(conSet.c_max)
# end
#
# function reset!(con::ConstraintVals{T,W,C,P}, opts::StaticALGamesSolverOptions{T}) where {T,W,C,P}
# 	λ = con.λ
# 	c = con.vals
# 	μ = con.μ
# 	for i in eachindex(con.inds)
# 		μ[i] = opts.penalty_initial * @SVector ones(T,P)
# 		c[i] *= 0.0
# 		λ[i] *= 0.0
# 	end
# end
