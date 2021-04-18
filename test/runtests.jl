using Test
using Algames
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
include("dynamics/double_integrator.jl")
include("dynamics/unicycle.jl")
include("dynamics/bicycle.jl")

# Struct
include("struct/problem_size.jl")
include("struct/statistics.jl")
include("struct/primal_dual_traj.jl")
include("struct/violations.jl")
include("struct/regularizer.jl")
include("struct/options.jl")

# Core
include("core/stamp.jl")
include("core/newton_core.jl")

# Constraints
include("constraints/state_bound_constraint.jl")
include("constraints/control_bound_constraint.jl")
include("constraints/wall_constraint.jl")
include("constraints/game_constraints.jl")
include("constraints/constraints_methods.jl")
include("constraints/constraint_derivatives.jl")

# Objective
include("objective/objective.jl")

# Problem
include("problem/problem.jl")
include("problem/local_quantities.jl")
include("problem/global_quantities.jl")
include("problem/solver_methods.jl")

# Equilibrium Subspace
include("active_set/active_set_stamp.jl")
include("active_set/active_set_core.jl")
include("active_set/active_set_methods.jl")

# Plots
include("plots/solver_plots.jl")
