export
    GameProblem,
    change_integration,
    integration,
    states,
    controls,
    initial_trajectory!,
    initial_states!,
    initial_controls!,
    max_violation,
    num_constraints,
    get_constraints,
    change_integration,
    rollout!

# """$(TYPEDEF) Trajectory Optimization Problem.
# Contains the full definition of a trajectory optimization problem, including:
# * dynamics model (`Model`)
# * objective (`Objective`)
# * constraints (`ConstraintSet`)
# * initial and final states
# * Primal variables (state and control trajectories)
# * Discretization information: knot points (`N`), time step (`dt`), and total time (`tf`)
#
# # Constructors:
# ```julia
# Problem(model, obj, constraints, x0, xf, Z, N, tf) # defaults to RK3 integration
# Problem{Q}(model, obj, constraints, x0, xf, Z, N, tf) where Q<:QuadratureRule
# Problem(model, obj, xf, tf; x0, constraints, N, X0, U0, dt, integration)
# Problem{Q}(prob::Problem)  # change integration
# ```
# where `Z` is a trajectory (Vector of `KnotPoint`s)
#
# # Arguments
# * `model`: Dynamics model. Can be either `Discrete` or `Continuous`
# * `obj`: Objective
# * `X0`: Initial state trajectory. If omitted it will be initialized with NaNs, to be later overwritten by the solver.
# * `U0`: Initial control trajectory. If omitted it will be initialized with zeros.
# * `x0`: Initial state. Defaults to zeros.
# * `xf`: Final state. Defaults to zeros.
# * `dt`: Time step
# * `tf`: Final time. Set to zero to specify a time penalized problem.
# * `N`: Number of knot points. Defaults to 51, unless specified by `dt` and `tf`.
# * `integration`: One of the defined integration types to discretize the continuous dynamics model.
# Both `X0` and `U0` can be either a `Matrix` or a `Vector{Vector}`, but must be the same.
# At least 2 of `dt`, `tf`, and `N` need to be specified (or just 1 of `dt` and `tf`).
# """
struct GameProblem{Q<:QuadratureRule,T<:AbstractFloat,O<:TO.AbstractObjective}
    model::AbstractGameModel
    obj::Vector{O}
    constraints::ConstraintSet{T}
    x0::SVector
    xf::SVector
    Z::Traj
    N::Int
    tf::T
    function GameProblem{Q}(
            model::AbstractGameModel,
            obj::Vector{O},
            constraints::ConstraintSet,
            x0::SVector,
            xf::SVector,
            Z::Traj,
            N::Int,
            tf::T) where {Q,T,O}

            n,m,pl,p = size(model)
            @assert length(x0) == length(xf) == n
            @assert length(Z) == N
        new{Q,T,O}(
            model,
            obj,
            constraints,
            x0,
            xf,
            Z,
            N,
            tf)
    end
end

"Use RK3 as default integration"
GameProblem(model, obj, constraints, x0, xf, Z, N, tf) =
    GameProblem{RK3}(model, obj, constraints, x0, xf, Z, N, tf)

function GameProblem(model::L, obj::Vector{O}, xf::AbstractVector, tf;
        constraints=ConstraintSet(length(obj[1])),
        x0=zero(xf), N::Int=length(obj[1]),
        X0=[x0*NaN for k = 1:N],
        U0=[@SVector zeros(size(model)[2]) for k = 1:N-1],
        dt=fill(tf/(N-1),N),
        integration=TO.DEFAULT_Q) where {L,O}
    n,m = size(model)
    if dt isa Real
        dt = fill(dt,N)
    end
    if X0 isa AbstractMatrix
        X0 = [X0[:,k] for k = 1:size(X0,2)]
    end
    if U0 isa AbstractMatrix
        U0 = [U0[:,k] for k = 1:size(U0,2)]
    end
    Z = Traj(X0,U0,dt)

    GameProblem{integration}(model, obj, constraints, SVector{n}(x0), SVector{n}(xf),
        Z, N, tf)
end



# "$(TYPEDSIGNATURES)
# Get number of states, controls, and knot points"
Base.size(prob::GameProblem) = size(prob.model)..., prob.N

# """```julia
# integration(::Problem)
# integration(::DynamicsConstraint)
# ```
# Get the integration rule"""
integration(prob::GameProblem{Q}) where Q = Q

# "```julia
# controls(::Problem)
# controls(::AbstractSolver)
# controls(::Traj)
# ```
# Get the control trajectory
# "
controls(prob::GameProblem) = controls(prob.Z)

# "```julia
# states(::Problem)
# states(::AbstractSolver)
# states(::Traj)
# ```
# Get the state trajectory
# "
states(prob::GameProblem) = states(prob.Z)

# "```julia
# initial_trajectory!(::Problem, Z)
# initial_trajectory!(::AbstractSolver, Z)
# ```
# Copy the trajectory "
function initial_trajectory!(prob::GameProblem, Z::Traj)
    for k = 1:prob.N
        prob.Z[k].z = Z[k].z
    end
end

# "```julia
# initial_states!(::Union{Problem,AbstractSolver}, X0::Vector{<:AbstractVector})
# initial_states!(::Union{Problem,AbstractSolver}, X0::AbstractMatrix)
# ```
# Copy the state trajectory "
function initial_states!(prob::GameProblem, X0::Vector{<:AbstractVector})
    set_states!(prob.Z, X0)
end

function initial_states!(prob::GameProblem, X0::AbstractMatrix)
    X0 = [X0[:,k] for k = 1:size(X0,2)]
    set_states!(prob.Z, X0)
end

# "```julia
# initial_controls!(::Union{Problem,AbstractSolver}, U0::Vector{<:AbstractVector})
# initial_controls!(::Union{Problem,AbstractSolver}, U0::AbstractMatrx)
# ```
# Copy the control trajectory "
function initial_controls!(prob::GameProblem, U0::Vector{<:AbstractVector})
    TO.set_controls!(prob.Z, U0)
end

function initial_controls!(prob::GameProblem, u0::AbstractVector{<:Real})
    U0 = [TO.copy(u0) for k = 1:prob.N]
    initial_controls!(prob, U0)
end

# "```julia
# cost(::Problem)
# cost(::AbstractSolver)
# ```
# Compute the cost for the current trajectory"
@inline TO.cost(prob::GameProblem) = cost(prob.obj, prob.Z)

"Copy the problem"
function TO.copy(prob::GameProblem{Q}) where Q
    GameProblem{Q}(prob.model, copy(prob.obj), copy(prob.constraints), prob.x0, prob.xf,
        copy(prob.Z), prob.N, prob.tf)
end


function TO.max_violation(prob::GameProblem, Z::Traj=prob.Z)
    conSet = get_constraints(prob)
    evaluate!(conSet, Z)
    max_violation!(conSet)
    return maximum(conSet.c_max)
end

TO.num_constraints(prob::GameProblem) = get_constraints(prob).p

@inline TO.get_constraints(prob::GameProblem) = prob.constraints


# "```julia
# change_integration(prob::Problem, Q<:QuadratureRule)
# ```
# Change dynamics integration for the problem"
TO.change_integration(prob::GameProblem, ::Type{Q}) where Q<:QuadratureRule =
    GameProblem{Q}(prob)

function GameProblem{Q}(p::Problem) where Q
    GameProblem{Q}(p.model, p.obj, p.constraints, p.x0, p.xf, p.Z, p.N, p.tf)
end

@inline rollout!(prob::GameProblem) = rollout!(prob.model, prob.Z, prob.x0)

function GameProblem(p::GameProblem; model=p.model, obj=p.obj, constraints=p.constraints,
    x0=p.x0, xf=p.xf, t0=p.t0, tf=p.tf)
    Problem(model, obj, constraints, x0, xf, p.Z, p.N, t0, tf)
end
