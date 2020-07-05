export
    DirectGamesStats,
    DirectGamesSolverOptions,
    DirectGamesSolver,
    reset!,
    converged,
    get_trajectory,
    get_objective,
    get_model,
    get_initial_state,
    cost


@with_kw mutable struct DirectGamesStats{T}
    iterations::Int = 0
    iterations_total::Int = 0
    iterations_inner::Vector{Int} = zeros(Int,0)
    cost::Vector{Vector{T}} = [zeros(T,0)]
    dJ::Vector{Vector{T}} = [zeros(T,0)]
    cmax::Vector{T} = zeros(T,0)
    optimality_merit::Vector{Vector{T}} = [zeros(T,0)]
    optimality_merit_inf::Vector{Vector{T}} = [zeros(T,0)]
    H_cond::Vector{Vector{T}} = [zeros(T,0)]
    α::Vector{Vector{T}} = [zeros(T,0)]
    dJ_zero_counter::Int = 0
    runtime::T = 0.
end

function reset!(stats::DirectGamesStats, L=0, l=0, p=0)
    stats.iterations = 0
    stats.iterations_total = 0
    stats.iterations_inner = zeros(Int,L)
    stats.cost = [zeros(p) for j = 1:L]
    stats.dJ = [zeros(p) for j = 1:L]
    stats.cmax = zeros(L)
    stats.optimality_merit = [zeros(l) for j = 1:L]
    stats.optimality_merit_inf = [zeros(l) for j = 1:L]
    stats.H_cond = [zeros(l) for j = 1:L]
    stats.α = [zeros(l) for j = 1:L]
    stats.dJ_zero_counter = 0
    stats.runtime = 0.
end

# """$(TYPEDEF)
# Solver options for the iterative LQGames (DirectGames) solver.
# $(FIELDS)
# """
@with_kw mutable struct DirectGamesSolverOptions{T} <: TO.AbstractSolverOptions{T}
    # Options

    "Print summary at each iteration."
    verbose::Bool=false

    "Live plotting."
    live_plotting::Symbol=:off # :state, :control

    "Compute H_ condition number."
    record_condition::Bool=false

    "type of game theoretic equilibrium, Nash or Stackelberg."
    eq_type::Symbol = :nash # :nash :stackelberg

    "type of information pattern, Memoryless Perfect State (MPS) leads to a feedback equilibrium, Open-Loop (OL) leads to an open-loop equilibrium."
    info_pattern::Symbol = :feedback # :feedback :open_loop

    "Primals regularization."
    η::T = 1e-10

    "Initial step size parameters for line search."
    α_init::T = 1.0

    "Initial step size parameters for line search."
    β::T = 0.1

    "Initial step size parameters for line search."
    τ::T = 1.0/2

    "dJ < ϵ, cost convergence criteria for unconstrained solve or to enter outerloop for constrained solve."
    cost_tolerance::T = 1.0e-4

    "|g_|_1 < ϵ, optimality constraint tolerance."
    optimality_constraint_tolerance::T = 1.0e-2

    "|g_|_∞ < ϵ, optimality constraint tolerance infinity norm."
    optimality_constraint_tolerance_inf::T = 1.0e+1

    "cmax < ϵ, constraint tolerance."
    constraint_tolerance::T = 1.0e-3

    "DirectGames iterations."
    iterations::Int = 10

    "DirectGames inner iterations."
    inner_iterations::Int = 20

    "restricts the total number of times a forward pass fails, resulting in regularization, before exiting."
    dJ_counter_limit::Int = 10

    "forward pass approximate line search lower bound, 0 < line_search_lower_bound < line_search_upper_bound."
    line_search_lower_bound::T = 1.0e-8

    "forward pass approximate line search upper bound, 0 < line_search_lower_bound < line_search_upper_bound < ∞."
    line_search_upper_bound::T = 10.0

    "maximum number of backtracking steps during forward pass line search."
    iterations_linesearch::Int = 10

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

    # Solver Numerical Limits
    "maximum cost value, if exceded solve will error."
    max_cost_value::T = 1.0e8

    "maximum state value, evaluated during rollout, if exceded solve will error."
    max_state_value::T = 1.0e8

    "maximum control value, evaluated during rollout, if exceded solve will error."
    max_control_value::T = 1.0e8

    "Minimum number of outer loops."
    min_iterations::Int = 0

    "Minimum number of inner steps taken in the inner loop before kicking out to the outer loop."
    min_steps_per_iteration::Int = 1

    "Allow the solver to break early out of the inner loop before making a step."
    break_inner_loop::Bool = true

    "Penalty scaling term"
    μ_penalty::T = 1.0

    "Time after which the solver kicks out and returns its current solution."
    timeout::T = Inf

    log_level::Base.CoreLogging.LogLevel = TO.InnerLoop
end

# """$(TYPEDEF)
# DirectGames is an unconstrained indirect method for trajectory optimization that parameterizes only the controls and enforces strict dynamics feasibility at every iteration by simulating forward the dynamics with an LQR feedback controller.
# The main algorithm consists of two parts:
# 1) a backward pass that uses Differential Dynamic Programming to compute recursively a quadratic approximation of the cost-to-go, along with linear feedback and feed-forward gain matrices, `K` and `d`, respectively, for an LQR tracking controller, and
# 2) a forward pass that uses the gains `K` and `d` to simulate forward the full nonlinear dynamics with feedback.
# """

# struct DirectGamesSolver{T,I<:QuadratureRule,L,O,n,n̄,m,L1,L2,D,F,E1,E2} <: ConstrainedSolver{T}
struct DirectGamesSolver{T,I<:QuadratureRule,L,O,n,m,L1,D,E1,nmi,nm,mi,l2n2mi,l2n2m,l2nmi,l2nm,} <: ConstrainedSolver{T}
    # Model + Objective
    model::L
    obj::Vector{O}

    spu::Vector{SVector{mi,Int}}

    xinds::Vector{SVector{n,Int}}
    uinds::Vector{SVector{m,Int}}
    xinds_p::Vector{Vector{SVector{n,Int}}}
    uinds_p::Vector{Vector{SVector{mi,Int}}}
    uinds_H::Vector{Vector{SVector{mi,Int}}}
    νinds::Vector{Vector{SVector{n,Int}}}
    νinds_p::Vector{SVector{n,Int}}
    sinds::StaticInds{n,m,nm,l2nm,l2n2m}
    stage_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{nmi,Int},SVector{nm,Int}},false}}}
    state_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{n,Int},SVector{n,Int}},false}}}
    control_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{mi,Int},SVector{m,Int}},false}}}
    coupled_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{l2n2mi,Int},SVector{l2n2m,Int}},false}}}
    dynamical_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{l2nmi,Int},SVector{l2nm,Int}},false}}}
    # state_view_g::Vector{Vector{SubArray{T,1,Vector{T},Tuple{SVector{n,Int}},false}}}

    # Problem info
    x0::SVector{n,T}
    xf::SVector{n,T}
    tf::T
    N::Int

    opts::DirectGamesSolverOptions{T}
    stats::DirectGamesStats{T}

    # Constraints
    constraints::ConstraintSet{T}
    penalty_constraints::ConstraintSet{T}
    dyn_constraints::ConstraintSet{T}

    # Primal Duals
    Z::Vector{KnotPoint{T,n,m,L1}}
    Z̄::Vector{KnotPoint{T,n,m,L1}}

    ∇F::Vector{D} # discrete dynamics jacobian (block) (n,n+m+1,N)
    C::Vector{E1}  # stage cost expansions trajectory

    ρ::Vector{T} # Regularization
    dρ::Vector{T} # Regularization rate of change
    η::Vector{T} # Primal regularization

    H_::SparseMatrixCSC{T,Int} # Pre-allocating
    g_::Vector{T} # Pre-allocating
    δY::Vector{T} # Pre-allocating

    ν::Vector{Vector{SVector{n,T}}}
    ν_::Vector{Vector{SVector{n,T}}}
    γ::Vector{Vector{SVector{n,T}}}

    logger::TO.SolverLogger

    function DirectGamesSolver{T,I}(model::L, obj::Vector{O},
            spu::Vector{SVector{mi,Int}},
            xinds::Vector{SVector{n,Int}},
            uinds::Vector{SVector{m,Int}},
            xinds_p::Vector{Vector{SVector{n,Int}}},
            uinds_p::Vector{Vector{SVector{mi,Int}}},
            uinds_H::Vector{Vector{SVector{mi,Int}}},
            νinds::Vector{Vector{SVector{n,Int}}},
            νinds_p::Vector{SVector{n,Int}},
            sinds::StaticInds{n,m,nm,l2nm,l2n2m},
            stage_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{nmi,Int},SVector{nm,Int}},false}}},
            state_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{n,Int},SVector{n,Int}},false}}},
            control_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{mi,Int},SVector{m,Int}},false}}},
            coupled_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{l2n2mi,Int},SVector{l2n2m,Int}},false}}},
            dynamical_view_H::Vector{Vector{SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{SVector{l2nmi,Int},SVector{l2nm,Int}},false}}},
            # state_view_g::Vector{Vector{SubArray{T,1,Vector{T},Tuple{SVector{n,Int}},false}}},
            x0,
            xf,
            tf,
            N,
            opts,
            stats,
            conSet::ConstraintSet{T},
            penalty_conSet::ConstraintSet{T},
            dyn_conSet::ConstraintSet{T},
            Z::Vector{KnotPoint{T,n,m,L1}},
            Z̄,
            ∇F::Vector{D},
            C::Vector{E1},
            ρ,
            dρ,
            η,
            H_::SparseMatrixCSC{T,Int},
            g_::Vector{T},
            δY::Vector{T},
            ν::Vector{Vector{SVector{n,T}}},
            ν_::Vector{Vector{SVector{n,T}}},
            γ::Vector{Vector{SVector{n,T}}},
            logger) where {T,I,L,O,n,m,L1,D,E1,nmi,nm,mi,l2n2mi,l2n2m,l2nmi,l2nm}
        new{T,I,L,O,n,m,L1,D,E1,nmi,nm,mi,l2n2mi,l2n2m,l2nmi,l2nm}(model, obj,
            spu,
            xinds,
            uinds,
            xinds_p,
            uinds_p,
            uinds_H,
            νinds,
            νinds_p,
            sinds,
            stage_view_H,
            state_view_H,
            control_view_H,
            coupled_view_H,
            dynamical_view_H,
            x0,
            xf,
            tf,
            N,
            opts,
            stats,
            conSet,
            penalty_conSet,
            dyn_conSet,
            Z,
            Z̄,
            ∇F,
            C,
            ρ,
            dρ,
            η,
            H_,
            g_,
            δY,
            ν,
            ν_,
            γ,
            logger)
    end
end

function DirectGamesSolver(prob::GameProblem{I,T}, opts=DirectGamesSolverOptions{T}()) where {I,T}
    n,m,pu,p = size(prob)
    N = prob.N
    # Init solver statistics
    stats = DirectGamesStats{T}()

    # Add dynamics constraints
    prob = copy(prob)
    conSet = get_constraints(prob)
    penalty_conSet = ConstraintSet(n,m,prob.N)
    dyn_conSet = ConstraintSet(n,m,prob.N)
    add_dynamics_constraints!(dyn_conSet,prob)

    # Init solver results
    n̄ = TO.state_diff_size(prob.model)
    spu = [SVector{length(pu[i])}(pu[i]) for i in 1:p]

    xinds = [-n+(k-1)*(m+n) .+ SVector{n}([j for j in 1:n]) for k=1:N] #okkk
    xinds_p = [[-n*p+(k-1)*(m+n*p)+(i-1)*n .+ SVector{n}([j for j in 1:n]) for k=1:N] for i=1:p] #okkk
    uinds = [(k-1)*(m+n) .+ SVector{m}([j for j in 1:m]) for k=1:N] #okkk
    uinds_p = [[(k-1)*(m+n*p) .+ spu[i] for k=1:N] for i=1:p] #okkk
    uinds_H = [[(k-1)*(m+n) .+ spu[i] for k=1:N] for i=1:p] #okkk
    off = (n+m)*(N-1)
    νinds = [[off+(k-1)*n*p+(i-1)*n .+ SVector{n}([j for j in 1:n]) for k=1:N-1] for i=1:p]
    off_p = (n*p+m)*(N-1)
    νinds_p = [off_p+(k-1)*n .+ SVector{n}([j for j in 1:n]) for k=1:N-1]

    sinds = StaticInds(
        SVector{0}(zeros(Int,0)),
        SVector{n}(1:n),
        SVector{m}(1:m),
        SVector{n+m}(1:n+m),
        SVector{2*n+m}(1:2*n+m),
        SVector{2*n+2*m}(1:2*n+2*m))

    x0 = SVector{n}(prob.x0)
    xf = SVector{n}(prob.xf)

    Z = prob.Z
    Z̄ = TO.copy(prob.Z)


    ∇F = [SMatrix{n,n+m+1}(zeros(T,n,n+m+1)) for k = 1:N-1]

    # Cost expansion of different sizes depending on the number of control for each player
    C = [TO.CostExpansion(n̄,length(pu[i]),N) for i=1:p]

    ρ = zeros(T,1)
    dρ = zeros(T,1)
    η = [opts.η]

    NN = (n*(p+1) + m)*(N-1)
    H_ = spzeros(T,NN,NN)
    g_ = zeros(T,NN)
    δY = zeros(T,NN)

    stage_view_H = [[view(H_,
        [xinds_p[i][k]; uinds_p[i][k]],
        [xinds[k]; uinds[k]]
        ) for k=2:N-1] for i=1:p]
    state_view_H = [[view(H_,
        xinds_p[i][k],
        xinds[k]
        ) for k=2:N] for i=1:p]
    control_view_H = [[view(H_,
        uinds_p[i][k],
        uinds[k]
        ) for k=1:N-1] for i=1:p]
    coupled_view_H = [[view(H_,
        [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]; uinds_p[i][k+1]],
        [xinds[k]; uinds[k]; xinds[k+1]; uinds[k+1]]
        ) for k=2:N-2] for i=1:p]
    dynamical_view_H = [[view(H_,
        [xinds_p[i][k]; uinds_p[i][k]; xinds_p[i][k+1]],
        [xinds[k]; uinds[k]; xinds[k+1]]
        ) for k=2:N-2] for i=1:p]

    # state_view_g = [[view(g_,xinds_p[i][k_1]) for k_1=2:N] for i=1:p]

    ν = [[SVector{n}(zeros(n)) for k=1:N-1] for i=1:p]
    ν_ = [[SVector{n}(zeros(n)) for k=1:N-1] for i=1:p]
    γ = [[SVector{n}(ones(n)) for k=1:N-1] for i=1:p]
    logger = TO.default_logger(opts.verbose)

    solver = DirectGamesSolver{T,I}(prob.model, prob.obj,
        spu,
        xinds,
        uinds,
        xinds_p,
        uinds_p,
        uinds_H,
        νinds,
        νinds_p,
        sinds,
        stage_view_H,
        state_view_H,
        control_view_H,
        coupled_view_H,
        dynamical_view_H,
        # state_view_g,
        x0,
        xf,
        prob.tf,
        N,
        opts,
        stats,
        conSet,
        penalty_conSet,
        dyn_conSet,
        Z,
        Z̄,
        ∇F,
        C,
        ρ,
        dρ,
        η,
        H_,
        g_,
        δY,
        ν,
        ν_,
        γ,
        logger)

    reset!(solver)
    return solver
end

AbstractSolver(prob::GameProblem, opts::DirectGamesSolverOptions) = DirectGamesSolver(prob, opts)

function reset!(solver::DirectGamesSolver{T}, reset_stats=true; reset_type::Symbol=:nominal) where T
    # The reset type can be either
    # - :nominal
    # - :full
    # - :mpc
    n,m,N = size(solver)
    if reset_stats
        n,m,pu,p = size(solver.model)
        reset!(solver.stats,
    		solver.opts.iterations,
    		solver.opts.inner_iterations,p)
    end
    TO.reset!(solver.constraints)
    TO.reset!(solver.penalty_constraints)
    nc = length(solver.penalty_constraints.constraints)
    set_penalty!(solver.penalty_constraints, solver.opts.μ_penalty[1])
    # if nc >= 1
    #     pen = solver.opts.μ_penalty[1] * ones(nc)
    #     set_penalty!(solver.penalty_constraints, pen)
    # end
    TO.reset!(solver.dyn_constraints)
    solver.ρ[1] = 0.0
    solver.dρ[1] = 0.0
    solver.η[1] = 0.0 #######
    if reset_type == :nominal || reset_type == :full
        for k = 1:N
            solver.Z̄[k].z = @SVector zeros(n+m)
            solver.Z[k].z = @SVector zeros(n+m)
        end
        for i = 1:p
            for k = 1:N-1
                solver.ν[i][k] = @SVector zeros(n)
            end
        end
    end
    if reset_type == :full
        solver.H_ .*= 0.0
        solver.g_ .*= 0.0
    end
    return nothing
end

function converged(solver::DirectGamesSolver{T,I}) where {T,I}
    iter = solver.stats.iterations
    inner_iter = solver.stats.iterations_inner[iter]
    out = (solver.stats.cmax[iter] <= solver.opts.constraint_tolerance) &&
        (solver.stats.optimality_merit[iter][inner_iter] <= solver.opts.optimality_constraint_tolerance) &&
        (solver.stats.optimality_merit_inf[iter][inner_iter] <= solver.opts.optimality_constraint_tolerance_inf)
    return out
end

Base.size(solver::DirectGamesSolver{T,I,L,O,n,m,L1,D,E1,nmi,nm,mi,l2n2mi,l2n2m,l2nmi,l2nm}) where {T,I,L,O,n,m,L1,D,E1,nmi,nm,mi,l2n2mi,l2n2m,l2nmi,l2nm} = n,m,solver.N
@inline TO.get_trajectory(solver::DirectGamesSolver) = solver.Z
@inline TO.get_objective(solver::DirectGamesSolver) = solver.obj
@inline TO.get_model(solver::DirectGamesSolver) = solver.model
@inline TO.get_initial_state(solver::DirectGamesSolver) = solver.x0

function TO.cost(solver::DirectGamesSolver, Z=solver.Z)
    n,m,pu,p = size(solver.model)
    n,m,N = size(solver)
    for i = 1:p
        cost!(solver.obj[i], Z, fill(pu[i],N))
    end
    return sum.(get_J.(solver.obj))
end

function DirectGamesSolver(solver_::DirectGamesSolver{T,I}, obj::Vector{O},
    x0::SVector{N,T}, xf::SVector{N,T}=solver_.xf, tf::T=solver_.tf) where {T,I,N,O}
    solver = DirectGamesSolver{T,I}(
        solver_.model,
        obj,
        solver_.spu,
        solver_.xinds,
        solver_.uinds,
        solver_.xinds_p,
        solver_.uinds_p,
        solver_.uinds_H,
        solver_.νinds,
        solver_.νinds_p,
        solver_.sinds,
        solver_.stage_view_H,
        solver_.state_view_H,
        solver_.control_view_H,
        solver_.coupled_view_H,
        solver_.dynamical_view_H,
        x0,
        xf,
        tf,
        solver_.N,
        solver_.opts,
        solver_.stats,
        solver_.constraints,
        solver_.penalty_constraints,
        solver_.dyn_constraints,
        solver_.Z,
        solver_.Z̄,
        solver_.∇F,
        solver_.C,
        solver_.ρ,
        solver_.dρ,
        solver_.η,
        solver_.H_,
        solver_.g_,
        solver_.δY,
        solver_.ν,
        solver_.ν_,
        solver_.γ,
        solver_.logger)
    return solver
end

function DirectGamesSolver(solver_::DirectGamesSolver{T,I}, conSet::ConstraintSet{T}=solver_.constraints) where {T,I}
    solver = DirectGamesSolver{T,I}(
        solver_.model,
        solver_.obj,
        solver_.spu,
        solver_.xinds,
        solver_.uinds,
        solver_.xinds_p,
        solver_.uinds_p,
        solver_.uinds_H,
        solver_.νinds,
        solver_.νinds_p,
        solver_.sinds,
        solver_.stage_view_H,
        solver_.state_view_H,
        solver_.control_view_H,
        solver_.coupled_view_H,
        solver_.dynamical_view_H,
        solver_.x0,
        solver_.xf,
        solver_.tf,
        solver_.N,
        solver_.opts,
        solver_.stats,
        conSet,
        solver_.penalty_constraints,
        solver_.dyn_constraints,
        solver_.Z,
        solver_.Z̄,
        solver_.∇F,
        solver_.C,
        solver_.ρ,
        solver_.dρ,
        solver_.η,
        solver_.H_,
        solver_.g_,
        solver_.δY,
        solver_.ν,
        solver_.ν_,
        solver_.γ,
        solver_.logger)
    return solver
end

function Base.copy(s::DirectGamesSolverOptions{T}) where T
    fnames = fieldnames(typeof(s))
    args = [getfield(s,fname) for fname in fnames]
    DirectGamesSolverOptions{T}(args...)
end
