module Algames

greet() = print("Hello World!")

using Altro
using BenchmarkTools
using ForwardDiff
using LinearAlgebra
using Parameters
using Printf
using Random
using RobotDynamics
using SparseArrays
using StaticArrays
using TrajectoryOptimization

# Dynamics
export
    AbstractGameModel,
    DoubleIntegratorGame,
    UnicycleGame,
    BicycleGame,
    dynamics

# Struct
export
    Options,
    ProblemSize,
    PrimalDualTraj,
    init_traj!,
    update_traj!,
    set_traj!,
    get_traj!,
    Δ_step,
    Statistics,
    record!,
    reset!,
    DynamicsViolation,
    dynamics_violation,
    ControlViolation,
    control_violation,
    StateViolation,
    state_violation,
    OptimalityViolation,
    optimality_violation,
    Regularizer,
    set!,
    mult!

# Constraints
export
    GameConstraintValues,
    ControlBoundConstraint,
    WallConstraint,
    add_collision_avoidance!,
    add_control_bound!,
    add_circle_constraint!,
    Wall,
    add_wall_constraint!,
    set_constraint_params!,
    reset!,
    reset_duals!,
    reset_penalties!,
    penalty_update!,
    dual_update!,
    dual_update,
    evaluate!,
    jacobian!,
    update_active_set!,
    constraint_jacobian_residual!,
    constraint_residual!

# Core
export
    NewtonCore,
    vertical_indices,
    horizontal_indices,
    idx,
    horizontal_idx,
    vertical_idx,
    residual_views,
    jacobian_views,
    dynamics_indices,
    Stamp,
    VStamp,
    HStamp,
    stampify,
    stampify!,
    valid

# Objective
export
    GameObjective,
    cost_gradient!,
    cost_hessian!,
    add_collision_cost!,
    CollisionCost

# Problem
export
    Penalty,
    GameProblem,
    residual!,
    residual_jacobian!,
    regularize_residual!, # need test
    regularize_residual_jacobian!, # need test
    add2sub,
    addI2sub,
    sparse_zero!,
    scn,
    dynamics_residual,
    ∇dynamics!,
    newton_solve!,
    inner_iteration,
    line_search

# Active Set
export
    NullSpace,
    add_matrix!,
    ActiveSetCore,
    complete_vertical_indices,
    complete_horizontal_indices,
    complete_residual_views,
    complete_jacobian_views,
    CStamp,
    active,
    active_vertical_mask!,
    active_horizontal_mask!,
    get_collision_conval,
    update_nullspace!


# Plots
export
    plot_traj!,
    plot_violation!

# Dynamics
include("dynamics/game_model.jl")
include("dynamics/double_integrator.jl")
include("dynamics/unicycle.jl")
include("dynamics/bicycle.jl")

# Struct
include("struct/problem_size.jl")

# Core
include("core/stamp.jl")
include("core/newton_core.jl")

# Struct
include("struct/primal_dual_traj.jl")
include("struct/regularizer.jl")
include("struct/options.jl")

# Constraints
include("constraints/control_bound_constraint.jl")
include("constraints/wall_constraint.jl")
include("constraints/game_constraints.jl")
include("constraints/constraints_methods.jl")

# Struct
include("struct/violations.jl")
include("struct/statistics.jl")

# Objective
include("objective/objective.jl")

# Problem
include("problem/problem.jl")
include("problem/local_quantities.jl")
include("problem/global_quantities.jl")
include("problem/solver_methods.jl")

# Constraints
include("constraints/constraint_derivatives.jl")

# Active Set
include("active_set/active_set_stamp.jl")
include("active_set/active_set_core.jl")
include("active_set/active_set_methods.jl")

# Plots
include("plots/solver_plots.jl")

end # module
