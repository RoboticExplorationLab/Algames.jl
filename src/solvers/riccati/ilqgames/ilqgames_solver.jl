export
    iLQGamesStats,
    iLQGamesSolverOptions,
    iLQGamesSolver,
    reset!,
    get_trajectory,
    get_objective,
    get_model,
    get_initial_state


@with_kw mutable struct iLQGamesStats{T}
    iterations::Int = 0
    cost::Vector{Vector{T}} = [zeros(p)]
    dJ::Vector{Vector{T}} = [zeros(p)]
    gradient::Vector{T} = [0.]
    dJ_zero_counter::Int = 0
end

function reset!(stats::iLQGamesStats, L=0, p=0)
    stats.iterations = 0
    stats.cost = [zeros(p) for l = 1:L]
    stats.dJ = [zeros(p) for l = 1:L]
    stats.gradient = zeros(N)
    stats.dJ_zero_counter = 0
end

# """$(TYPEDEF)
# Solver options for the iterative LQGames (iLQGames) solver.
# $(FIELDS)
# """
@with_kw mutable struct iLQGamesSolverOptions{T} <: TO.AbstractSolverOptions{T}
    # Options

    "Print summary at each iteration."
    verbose::Bool=false

    "Live plotting."
    live_plotting::Symbol=:off # :state, :control

    "type of game theoretic equilibrium, Nash or Stackelberg."
    eq_type::Symbol = :nash # :nash :stackelberg

    "type of information pattern, Memoryless Perfect State (MPS) leads to a feedback equilibrium, Open-Loop (OL) leads to an open-loop equilibrium."
    info_pattern::Symbol = :feedback # :feedback :open_loop

    "dJ < ϵ, cost convergence criteria for unconstrained solve or to enter outerloop for constrained solve."
    cost_tolerance::T = 1.0e-4

    "gradient type: :todorov, :feedforward."
    gradient_type::Symbol = :todorov

    "gradient_norm < ϵ, gradient norm convergence criteria."
    gradient_norm_tolerance::T = 1.0e-5

    "iLQGames iterations."
    iterations::Int = 300

    "restricts the total number of times a forward pass fails, resulting in regularization, before exiting."
    dJ_counter_limit::Int = 10

    "use square root method backward pass for numerical conditioning."
    square_root::Bool = false

    "forward pass approximate line search lower bound, 0 < line_search_lower_bound < line_search_upper_bound."
    line_search_lower_bound::T = 1.0e-8

    "forward pass approximate line search upper bound, 0 < line_search_lower_bound < line_search_upper_bound < ∞."
    line_search_upper_bound::T = 10.0

    "maximum number of backtracking steps during forward pass line search."
    iterations_linesearch::Int = 20

    # Regularization
    "initial regularization."
    bp_reg_initial::T = 0.0

    "regularization scaling factor."
    bp_reg_increase_factor::T = 1.6

    "maximum regularization value."
    bp_reg_max::T = 1.0e8

    "minimum regularization value."
    bp_reg_min::T = 1.0e-8

    "type of regularization- control: () + ρI, state: (S + ρI); see Synthesis and Stabilization of Complex Behaviors through Online Trajectory Optimization."
    bp_reg_type::Symbol = :control

    "additive regularization when forward pass reaches max iterations."
    bp_reg_fp::T = 10.0

    # square root backward pass options:
    "type of matrix inversion for bp sqrt step."
    bp_sqrt_inv_type::Symbol = :pseudo

    "initial regularization for square root method."
    bp_reg_sqrt_initial::T = 1.0e-6

    "regularization scaling factor for square root method."
    bp_reg_sqrt_increase_factor::T = 10.0

    # Solver Numerical Limits
    "maximum cost value, if exceded solve will error."
    max_cost_value::T = 1.0e8

    "maximum state value, evaluated during rollout, if exceded solve will error."
    max_state_value::T = 1.0e8

    "maximum control value, evaluated during rollout, if exceded solve will error."
    max_control_value::T = 1.0e8

    log_level::Base.CoreLogging.LogLevel = TO.InnerLoop
end


# """$(TYPEDEF)
# iLQGames is an unconstrained indirect method for trajectory optimization that parameterizes only the controls and enforces strict dynamics feasibility at every iteration by simulating forward the dynamics with an LQR feedback controller.
# The main algorithm consists of two parts:
# 1) a backward pass that uses Differential Dynamic Programming to compute recursively a quadratic approximation of the cost-to-go, along with linear feedback and feed-forward gain matrices, `K` and `d`, respectively, for an LQR tracking controller, and
# 2) a forward pass that uses the gains `K` and `d` to simulate forward the full nonlinear dynamics with feedback.
# """

struct iLQGamesSolver{T,I<:QuadratureRule,L,O,n,n̄,m,L1,L2,D,F,E1,E2} <: UnconstrainedSolver{T}
    # Model + Objective
    model::L
    obj::Vector{O}

    # Problem info
    x0::SVector{n,T}
    xf::SVector{n,T}
    tf::T
    N::Int

    opts::iLQGamesSolverOptions{T}
    stats::iLQGamesStats{T}

    # Primal Duals
    Z::Vector{KnotPoint{T,n,m,L1}}
    Z̄::Vector{KnotPoint{T,n,m,L1}}

    # Data variables
    K::Vector{SMatrix{m,n̄,T,L2}}  # State feedback gains (m,n,N-1)
    d::Vector{SVector{m,T}} # Feedforward gains (m,N-1)

    ∇F::Vector{D} # discrete dynamics jacobian (block) (n,n+m+1,N)
    G::Vector{F}  # state difference jacobian (n̄, n)

    S::Vector{E1}  # Optimal cost-to-go expansion trajectory
    Q::Vector{E2}  # cost-to-go expansion trajectory
    C::Vector{E1}  # stage cost expansions trajectory

    ρ::Vector{T} # Regularization
    dρ::Vector{T} # Regularization rate of change

    grad::Vector{T} # Gradient

    N_::Array{T,2} # Pre-allocating for backward pass
    M_::Array{T,2} # Pre-allocating for backward pass
    m_::Vector{T} # Pre-allocating for backward pass

    logger::TO.SolverLogger

    function iLQGamesSolver{T,I}(
        model::L,
        obj::Vector{O},
        x0,
        xf,
        tf,
        N,
        opts,
        stats,
        Z::Vector{KnotPoint{T,n,m,L1}},
        Z̄,
        K::Vector{SMatrix{m,n̄,T,L2}},
        d,
        ∇F::Vector{D},
        G::Vector{F},
        S::Vector{E1},
        Q::Vector{E2},
        C::Vector{E1},
        ρ,
        dρ,
        grad,
        N_::Array{T,2},
        M_::Array{T,2},
        m_::Vector{T},
        logger) where {T,I,L,O,n,n̄,m,L1,L2,D,F,E1,E2}
        new{T,I,L,O,n,n̄,m,L1,L2,D,F,E1,E2}(
            model,
            obj,
            x0,
            xf,
            tf,
            N,
            opts,
            stats,
            Z,
            Z̄,
            K,
            d,
            ∇F,
            G,
            S,
            Q,
            C,
            ρ,
            dρ,
            grad,
            N_,
            M_,
            m_,
            logger)
    end
end

function iLQGamesSolver(prob::GameProblem{I,T}, opts=iLQGamesSolverOptions()) where {I,T}

    # Init solver statistics
    stats = iLQGamesStats{T}() # = Dict{Symbol,Any}(:timer=>TimerOutput())

    # Init solver results
    n,m,pl,p,N = size(prob)
    n̄ = TO.state_diff_size(prob.model)

    x0 = SVector{n}(prob.x0)
    xf = SVector{n}(prob.xf)

    Z = prob.Z
    Z̄ = TO.copy(prob.Z)

    K  = [@SMatrix zeros(T,m,n̄) for k = 1:N-1]
    d  = [@SVector zeros(T,m)   for k = 1:N-1]

    ∇F = [@SMatrix zeros(T,n,n+m+1) for k = 1:N-1]
    G = [TO.state_diff_jacobian(prob.model, x0) for k = 1:N]

    S = [TO.CostExpansion(n̄,length(pl[i]),N) for i=1:p]
    # Cost expansion of different sizes depending on the number of control for each player
    Q = [TO.CostExpansion(n,m,N) for i=1:p]
    C = [TO.CostExpansion(n̄,length(pl[i]),N) for i=1:p]

    ρ = zeros(T,1)
    dρ = zeros(T,1)

    grad = zeros(T,N-1)

    N_ = zeros(T,m,m)
    M_ = zeros(T,m,n)
    m_ = zeros(T,m)

    logger = TO.default_logger(opts.verbose)

    solver = iLQGamesSolver{T,I}(prob.model, prob.obj, x0, xf, prob.tf, N, opts, stats,
        Z, Z̄, K, d, ∇F, G, S, Q, C, ρ, dρ, grad, N_, M_, m_, logger)

    reset!(solver)
    return solver
end

AbstractSolver(prob::GameProblem, opts::iLQGamesSolverOptions) = iLQGamesSolver(prob, opts)

# function reset!(solver::iLQGamesSolver{T}, reset_stats=true) where T
#     if reset_stats
#         reset!(solver.stats, solver.opts.iterations)
#     end
#     n,m,N = size(solver)
#     for k = 1:N
#         solver.Z̄[k].z = @SVector zeros(n+m)
#         solver.Z[k].z = @SVector zeros(n+m)
#     end
    # TO.set_state!(solver.Z̄[1], solver.x0)
    # TO.set_state!(solver.Z[1], solver.x0)
#     solver.ρ[1] = 0.0
#     solver.dρ[1] = 0.0
#     return nothing
# end


function reset!(solver::iLQGamesSolver{T}, reset_stats=true; reset_type::Symbol=:nominal) where T
    # The reset type can be either
    # - :nominal
    # - :full
    # - :mpc
    n,m,N = size(solver)
    n,m,pl,p = size(solver.model)
    if reset_stats
        reset!(solver.stats, solver.opts.iterations, p)
    end
    solver.ρ[1] = 0.0
    solver.dρ[1] = 0.0
    if reset_type == :nominal || reset_type == :full
        for k = 1:N
            solver.Z̄[k].z = SVector{n+m}(zeros(n+m))
            solver.Z[k].z = SVector{n+m}(zeros(n+m))
        end
    end
    # TO.set_state!(solver.Z̄[1], solver.x0)
    # TO.set_state!(solver.Z[1], solver.x0)
    return nothing
end



Base.size(solver::iLQGamesSolver{T,I,L,O,n,m}) where {T,I,L,O,n,m} = solver.model.n,solver.model.m,solver.N
@inline TO.get_trajectory(solver::iLQGamesSolver) = solver.Z
@inline TO.get_objective(solver::iLQGamesSolver) = solver.obj
@inline TO.get_model(solver::iLQGamesSolver) = solver.model
@inline TO.get_initial_state(solver::iLQGamesSolver) = solver.x0

function TO.cost(solver::iLQGamesSolver, Z=solver.Z)
    n,m,pl,p = size(solver.model)
    n,m,N = size(solver)
    for i = 1:p
        cost!(solver.obj[i], Z, fill(pl[i],N))
    end
    return sum.(get_J.(solver.obj))
end
