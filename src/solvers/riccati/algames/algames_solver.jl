export
    ALGamesStats,
    ALGamesSolverOptions,
    ALGamesSolver,
    reset!,
    set_verbosity!,
    cost,
    get_trajectory,
    get_objective,
    get_model,
    get_initial_state,
    get_constraints,
    cost!,
    cost_expansion!

@with_kw mutable struct ALGamesStats{T}
    iterations::Int = 0
    iterations_total::Int = 0
    iterations_inner::Vector{Int} = zeros(Int,0)
    cost::Vector{Vector{T}} = [zeros(p)]
    c_max::Vector{T} = zeros(0)
    penalty_max::Vector{T} = zeros(0)
end

function reset!(stats::ALGamesStats, L=0, p=0)
    stats.iterations = 0
    stats.iterations_total = 0
    stats.iterations_inner = zeros(Int,L)
    stats.cost = [zeros(L)*NaN for i=1:p]
    stats.c_max = zeros(L)*NaN
    stats.penalty_max = zeros(L)*NaN
end


# """$(TYPEDEF)
# Solver options for the augmented Lagrangian solver.
# $(FIELDS)
# """
@with_kw mutable struct ALGamesSolverOptions{T} <: TO.AbstractSolverOptions{T}
    "Print summary at each iteration."
    verbose::Bool=false

    "unconstrained solver options."
    opts_uncon::TO.AbstractSolverOptions{T} = iLQGamesSolverOptions{Float64}()

    "dJ < ϵ, cost convergence criteria for unconstrained solve or to enter outerloop for constrained solve."
    cost_tolerance::T = 1.0e-4

    "dJ < ϵ_int, intermediate cost convergence criteria to enter outerloop of constrained solve."
    cost_tolerance_intermediate::T = 1.0e-3

    "gradient_norm < ϵ, gradient norm convergence criteria."
    gradient_norm_tolerance::T = 1.0e-5

    "gradient_norm_int < ϵ, gradient norm intermediate convergence criteria."
    gradient_norm_tolerance_intermediate::T = 1.0e-5

    "max(constraint) < ϵ, constraint convergence criteria."
    constraint_tolerance::T = 1.0e-3

    "max(constraint) < ϵ_int, intermediate constraint convergence criteria."
    constraint_tolerance_intermediate::T = 1.0e-3

    "maximum outerloop updates."
    iterations::Int = 30

    "global maximum Lagrange multiplier. If NaN, use value from constraint"
    dual_max::T = NaN

    "global maximum penalty term. If NaN, use value from constraint"
    penalty_max::T = NaN

    "global initial penalty term. If NaN, use value from constraint"
    penalty_initial::T = NaN

    "global penalty update multiplier; penalty_scaling > 1. If NaN, use value from constraint"
    penalty_scaling::T = NaN

    "penalty update multiplier when μ should not be update, typically 1.0 (or 1.0 + ϵ)."
    penalty_scaling_no::T = 1.0

    "ratio of current constraint to previous constraint violation; 0 < constraint_decrease_ratio < 1."
    constraint_decrease_ratio::T = 0.25

    "type of outer loop update (default, feedback)."
    outer_loop_update_type::Symbol = :default

    "numerical tolerance for constraint violation."
    active_constraint_tolerance::T = 0.0

    "terminal solve when maximum penalty is reached."
    kickout_max_penalty::Bool = false

    log_level::Base.CoreLogging.LogLevel = TO.OuterLoop
end

function reset!(conSet::ConstraintSet{T}, opts::ALGamesSolverOptions{T}) where T
    TO.reset!(conSet)
    if !isnan(opts.dual_max)
        for con in conSet.constraints
            params = get_params(con)::ConstraintParams{T}
            params.λ_max = opts.dual_max
        end
    end
    if !isnan(opts.penalty_max)
        for con in conSet.constraints
            params = get_params(con)::ConstraintParams{T}
            params.μ_max = opts.penalty_max
        end
    end
    if !isnan(opts.penalty_initial)
        for con in conSet.constraints
            params = get_params(con)::ConstraintParams{T}
            params.μ0 = opts.penalty_initial
        end
    end
    if !isnan(opts.penalty_scaling)
        for con in conSet.constraints
            params = get_params(con)::ConstraintParams{T}
            params.ϕ = opts.penalty_scaling
        end
    end
end

function set_verbosity!(opts::ALGamesSolverOptions)
    log_level = opts.log_level
    if opts.verbose
        set_logger()
        Logging.disable_logging(LogLevel(log_level.level-1))
        logger = global_logger()
        if opts.opts_uncon.verbose
            freq = 1
        else
            freq = 5
        end
        logger.leveldata[log_level].freq = freq
    else
        Logging.disable_logging(log_level)
    end
end


# @doc raw""" ```julia
# struct AugmentedLagrangianSolver <: TrajectoryOptimization.AbstractSolver{T}
# ```
# Augmented Lagrangian (AL) is a standard tool for constrained optimization. For a trajectory optimization problem of the form:
# ```math
# \begin{aligned}
#   \min_{x_{0:N},u_{0:N-1}} \quad & \ell_f(x_N) + \sum_{k=0}^{N-1} \ell_k(x_k, u_k, dt) \\
#   \textrm{s.t.}            \quad & x_{k+1} = f(x_k, u_k), \\
#                                  & g_k(x_k,u_k) \leq 0, \\
#                                  & h_k(x_k,u_k) = 0.
# \end{aligned}
# ```
# AL methods form the following augmented Lagrangian function:
# ```math
# \begin{aligned}
#     \ell_f(x_N) + &λ_N^T c_N(x_N) + c_N(x_N)^T I_{\mu_N} c_N(x_N) \\
#            & + \sum_{k=0}^{N-1} \ell_k(x_k,u_k,dt) + λ_k^T c_k(x_k,u_k) + c_k(x_k,u_k)^T I_{\mu_k} c_k(x_k,u_k)
# \end{aligned}
# ```
# This function is then minimized with respect to the primal variables using any unconstrained minimization solver (e.g. iLQR).
#     After a local minima is found, the AL method updates the Lagrange multipliers λ and the penalty terms μ and repeats the unconstrained minimization.
#     AL methods have superlinear convergence as long as the penalty term μ is updated each iteration.
# """
struct ALGamesSolver{T,S<:TO.AbstractSolver} <: TO.ConstrainedSolver{T}
    opts::ALGamesSolverOptions{T}
    stats::ALGamesStats{T}
    stats_uncon::Vector{STATS} where STATS
    solver_uncon::S
end

AbstractSolver(prob::GameProblem{Q,T},
    opts::ALGamesSolverOptions{T}=ALGamesSolverOptions{T}()) where {Q,T} =
    ALGamesSolver(prob,opts)

# """$(TYPEDSIGNATURES)
# Form an augmented Lagrangian cost function from a Problem and AugmentedLagrangianSolver.
#     Does not allocate new memory for the internal arrays, but points to the arrays in the solver.
# """
function ALGamesSolver(prob::GameProblem{Q,T}, opts::ALGamesSolverOptions=ALGamesSolverOptions{T}()) where {Q,T}
    # Init solver statistics
    stats = ALGamesStats()
    stats_uncon = Vector{iLQGamesSolverOptions{T}}()

    # Convert problem to AL problem
    n,m,pl,p = size(prob.model)
    N = prob.N
    alobj = [ALObjective(prob.obj[i], prob.constraints) for i=1:p]
    rollout!(prob)

    prob_al = GameProblem(prob.model, alobj, ConstraintSet(n,m,N),
        prob.x0, prob.xf, prob.Z, prob.N, prob.tf)

    solver_uncon = AbstractSolver(prob_al, opts.opts_uncon)

    solver = ALGamesSolver(opts,stats,stats_uncon,solver_uncon)
    reset!(solver)
    return solver
end

function reset!(solver::ALGamesSolver{T}) where T
    reset!(solver.stats, solver.opts.iterations, solver.solver_uncon.model.p)
    reset!(solver.solver_uncon)
    for con in TO.get_constraints(solver)
        reset!(con, solver.opts)
    end
end


Base.size(solver::ALGamesSolver) = size(solver.solver_uncon)
@inline TO.cost(solver::ALGamesSolver) = TO.cost(solver.solver_uncon)
@inline TO.get_trajectory(solver::ALGamesSolver) = TO.get_trajectory(solver.solver_uncon)
@inline TO.get_objective(solver::ALGamesSolver) = TO.get_objective(solver.solver_uncon)
@inline TO.get_model(solver::ALGamesSolver) = TO.get_model(solver.solver_uncon)
@inline TO.get_initial_state(solver::ALGamesSolver) = TO.get_initial_state(solver.solver_uncon)



function TO.get_constraints(solver::ALGamesSolver{T}) where T
    @show typeof(solver.solver_uncon)
    @show typeof(TO.get_objective(solver.solver_uncon))
    @show typeof(TO.get_objective(solver))
    obj = TO.get_objective(solver)#::ALObjective{T}
    # obj = TO.get_objective(solver)::ALObjective{T}
    [obj_i.constraints for obj_i in obj]
end




struct ALObjective{T,O<:Objective} <: TO.AbstractObjective
    obj::O
    constraints::ConstraintSet{T}
end

get_J(obj::ALObjective) = obj.obj.J
Base.length(obj::ALObjective) = length(obj.obj)

# TrajectoryOptimization.num_constraints(prob::Problem{Q,T,<:ALObjective}) where {T,Q} = prob.obj.constraints.p

function Base.copy(obj::ALObjective)
    ALObjective(obj.obj, ConstraintSet(copy(obj.constraints.constraints), length(obj.obj)))
end

# function cost!(obj::ALObjective, Z::Traj)
#     # Calculate unconstrained cost
#     cost!(obj.obj, Z)
#
#     # Calculate constrained cost
#     evaluate!(obj.constraints, Z)
#     update_active_set!(obj.constraints, Z, Val(0.0))
#     for con in obj.constraints.constraints
#         cost!(obj.obj.J, con, Z)
#     end
# end


function cost!(obj::ALObjective, Z::Traj, list_ind::Vector{Vector{Int}})
    # Calculate unconstrained cost
    cost!(obj.obj, Z, list_inds)

    # Calculate constrained cost
    evaluate!(obj.constraints, Z)
    @warn "need to have a set of constraint for each player so that
        they are not mixed when evaluating the indivudial costs"
    update_active_set!(obj.constraints, Z, Val(0.0))
    for con in obj.constraints.constraints
        cost!(obj.obj.J, con, Z)
    end
end

function cost_expansion!(E, obj::ALObjective, Z::Traj)
    # Update constraint jacobians
    jacobian!(obj.constraints, Z)

    ix, iu = Z[1]._x, Z[1]._u

    # Calculate expansion of original objective
    cost_expansion!(E, obj.obj, Z)

    # Add in expansion of constraints
    for con in obj.constraints.constraints
        cost_expansion(E, con, Z)
    end
end
