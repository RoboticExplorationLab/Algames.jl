"""
    GameSolver
Primary module for setting up and solving dynamic games problems.
"""
module GameSolver

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

export
	MonteCarloSampler,
	MonteCarloSamplerStats,
	MonteCarloSamplerOptions,
	reset!,
	get_trajectory,
	get_objective,
	get_model,
	get_initial_state

export
	monte_carlo_sampling

export
    Scenario,
    IntersectionScenario,
    StraightScenario,
    TIntersectionScenario,
    MergingScenario,
	add_rounded_boundary_constraints,
	add_scenario_constraints

export
	scope,
	animation,
	build_actors,
	scene_animation,
	actor_transformations

export
	plan_animation,
	animation,
	build_trajectory,
	still_scene_animation,
	still_animation

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
    MPCGamesSolverOptions,
    MPCGamesSolver

export
    ALGamesSolver,
    ALGamesSolverOptions,
    get_constraints

# iLQGames Solver
export
    iLQGamesSolverOptions,
    iLQGamesSolver

# iLQGamesPenalty Solver
export
    PenaltyiLQGamesSolverOptions,
    PenaltyiLQGamesSolver

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
    change_integration

export
	rel_zinds,
	zinds

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
    test_allocation,
    timing_solve


include("src/solvers/game_model.jl")
include("src/solvers/game_problem.jl")
include("src/solvers/cost_helpers.jl")

include("src/solvers/direct/direct_helpers.jl")

include("src/solvers/direct/direct_solver.jl")
include("src/solvers/direct/direct_methods.jl")
include("src/solvers/direct/direct_core.jl")
include("src/solvers/direct/newton_gradient.jl")
include("src/solvers/direct/newton_hessian.jl")
include("src/solvers/inds_helpers.jl")


include("src/solvers/riccati/algames/algames_solver.jl")
include("src/solvers/riccati/algames/algames_methods.jl")
include("src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("src/solvers/riccati/ilqgames/ilqgames_methods.jl")
include("src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("src/sampler/monte_carlo_sampler.jl")
include("src/sampler/monte_carlo_methods.jl")

include("src/scenarios/scenario.jl")
include("src/scenarios/examples/merging.jl")
include("src/scenarios/examples/straight.jl")
include("src/scenarios/examples/t_intersection.jl")

include("src/solvers/MPC/mpc_solver.jl")
include("src/solvers/MPC/mpc_methods.jl")

include("src/scenarios/scenario_mpc_visualization.jl")
include("src/scenarios/scenario_visualization.jl")

include("src/utils/constraints.jl")
include("src/utils/monte_carlo_visualization_latex.jl")
include("src/utils/monte_carlo_visualization.jl")
include("src/utils/plot_visualization.jl")
include("src/utils/tests.jl")
include("src/utils/timing.jl")

end
