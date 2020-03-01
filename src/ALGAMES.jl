"""
    ALGAMES
Primary module for setting up and solving dynamic games problems.
"""
module ALGAMES

using BenchmarkTools
using Blink
using Colors: RGBA, RGB
using CoordinateTransformations
using Dates
using FileIO
using GeometryTypes
using JLD2
using LinearAlgebra
using Logging
using MeshCat
using MeshIO
using Parameters
using PartedArrays
using PGFPlotsX
using Plots
using Random
using SparseArrays
using StaticArrays
using Statistics
using StatsBase
using Test
using TrajectoryOptimization
const TO = TrajectoryOptimization

using TrajectoryOptimization.Dynamics
using TrajectoryOptimization.Problems

# Sampler
export
	MonteCarloSampler,
	MonteCarloSamplerStats,
	MonteCarloSamplerOptions,
	reset!,
	record_sample,
	get_trajectory,
	get_objective,
	get_model,
	get_initial_state

export
	monte_carlo_sampling

# Scenario
export
    Scenario,
    IntersectionScenario,
    StraightScenario,
    TIntersectionScenario,
    MergingScenario,
	add_rounded_boundary_constraints,
	add_scenario_constraints

export
	solver_scope,
	animation,
	build_actors,
	build_trajectory,
	scene_animation,
	actor_transformations

export
    plot_velocity,
    latex_plot_velocity

export
	build_scenario,
	add_scenario_constraints,
	add_merging_constraints,
	add_straight_constraints,
	add_t_intersection_constraints,
	t_intersection_landmarks,
	add_rounded_boundary_constraints

# DirectGames Solver
export
	dual_ascent!,
	penalty_update!,
	inner_step!,
	line_search!,
	primal_dual_copy_update!,
	primal_dual_update!,
	indiv_cost
export
	add_dynamics_constraints!,
	StaticInds
export
	solve!,
	step!,
	record_iteration!,
	record_inner_iteration!,
	evaluate_convergence,
	evaluate_inner_convergence,
	rollout!,
	regularization_update!,
	regularize_primals!
export
    DirectGamesStats,
    DirectGamesSolverOptions,
    DirectGamesSolver,
    reset!,
    get_trajectory,
    get_objective,
    get_model,
    get_initial_state,
    cost
export
	update_g_!,
	update_g_con!
export
	update_H_!,
	update_H_con!,
	set_view_H!

# MPC Solver
export
	solve!,
	update_traj!,
	step!,
	need_resolve,
	generate_new_x0,
	mpc_propagate_dynamics,
	record_iteration!,
	evaluate_convergence,
	resample!

export
	MPCGamesStats,
	MPCGamesSolverOptions,
	MPCGamesSolver,
	reset!,
	get_trajectory,
	get_objective,
	get_model,
	get_initial_state

# ALGAMES
# export
# 	solve!,
# 	step!,
# 	record_iteration!,
# 	set_tolerances!,
# 	evaluate_convergence,
# 	dual_update!,
# 	penalty_update!

# export
#     ALGamesStats,
#     ALGamesSolverOptions,
#     ALGamesSolver,
#     reset!,
#     set_verbosity!,
#     cost,
#     get_trajectory,
#     get_objective,
#     get_model,
#     get_initial_state,
#     get_constraints,
#     cost!,
#     cost_expansion!

# iLQGames Solver
export
    solve!,
    step!,
    record_iteration!,
    gradient_todorov!,
    evaluate_convergence,
    nash_feedback_backwardpass!,
    nash_open_loop_backwardpass!,
    forwardpass!,
    rollout!,
    regularization_update!

export
    iLQGamesStats,
    iLQGamesSolverOptions,
    iLQGamesSolver,
    reset!,
    get_trajectory,
    get_objective,
    get_model,
    get_initial_state

# iLQGamesPenalty Solver
export
	solve!,
	set_penalty!,
	step!,
	penalty_expansion!,
	record_iteration!,
	gradient_todorov!,
	evaluate_convergence,
	nash_feedback_backwardpass!,
	nash_open_loop_backwardpass!,
	forwardpass!,
	rollout!,
	regularization_update!

export
    PenaltyiLQGamesSolverOptions,
    PenaltyiLQGamesStats,
    PenaltyiLQGamesSolver,
    reset!,
    get_trajectory,
    get_objective,
    get_constraints,
    get_model,
    get_initial_state

export
    stage_cost,
    cost!,
    cost_expansion,
    cost_gradient,
    cost_hessian,
    cost_gradient!,
    cost_hessian!

export
	AbstractGameModel

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

export
	rel_zinds,
	zinds

# Utils
export
	CollisionConstraint,
	BoundaryConstraint,
	ExpCircleConstraint,
	add_collision_avoidance,
	add_leader_constraints,
	add_circle_boundary,
	state_dim,
	evaluate

export
	visualize_latex_cmax,
	visualize_latex_H_cond,
	visualize_latex_solve_time,
	visualize_latex_iterations_total,
	visualize_latex_optimality_merit,
	visualize_latex_sampler

export
    visualize_cmax,
    visualize_H_cond,
    visualize_solve_time,
    visualize_iterations_total,
    visualize_optimality_merit

export
    visualize_trajectory_car,
    visualize_control,
    visualize_state,
    visualize_collision_avoidance,
    visualize_circle_collision,
    visualize_boundary_collision,
    visualize_dynamics,
    visualize_optimality_merit,
    visualize_H_cond,
    visualize_α,
    visualize_cmax,
    circle_shape,
    rectangle_shape

export
    test_allocation

export
    timing_solve



include("solvers/game_model.jl")
include("solvers/game_problem.jl")
include("solvers/cost_helpers.jl")

include("solvers/direct/direct_helpers.jl")

include("solvers/direct/direct_solver.jl")
include("solvers/direct/direct_methods.jl")
include("solvers/direct/direct_core.jl")
include("solvers/direct/newton_gradient.jl")
include("solvers/direct/newton_hessian.jl")
include("solvers/inds_helpers.jl")


include("solvers/riccati/algames/algames_solver.jl")
include("solvers/riccati/algames/algames_methods.jl")
include("solvers/riccati/ilqgames/ilqgames_solver.jl")
include("solvers/riccati/ilqgames/ilqgames_methods.jl")
include("solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("sampler/monte_carlo_sampler.jl")
include("sampler/monte_carlo_methods.jl")

include("scenarios/scenario.jl")
include("scenarios/examples/merging.jl")
include("scenarios/examples/straight.jl")
include("scenarios/examples/t_intersection.jl")

include("solvers/MPC/mpc_solver.jl")
include("solvers/MPC/mpc_methods.jl")

include("scenarios/scenario_visualization.jl")
include("scenarios/adaptive_plot.jl")

include("utils/constraints.jl")
include("utils/monte_carlo_visualization_latex.jl")
include("utils/monte_carlo_visualization.jl")
include("utils/plot_visualization.jl")
include("utils/tests.jl")
include("utils/timing.jl")

end
