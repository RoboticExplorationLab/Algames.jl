
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

include("../../src/solvers/game_model.jl")
include("../../src/solvers/game_problem.jl")
include("../../src/solvers/cost_helpers.jl")

include("../../src/solvers/direct/direct_helpers.jl")

include("../../src/solvers/direct/direct_solver.jl")
include("../../src/solvers/direct/direct_methods.jl")
include("../../src/solvers/direct/direct_core.jl")
include("../../src/solvers/direct/newton_gradient.jl")
include("../../src/solvers/direct/newton_hessian.jl")
include("../../src/solvers/inds_helpers.jl")


include("../../src/solvers/riccati/algames/algames_solver.jl")
include("../../src/solvers/riccati/algames/algames_methods.jl")
include("../../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("../../src/solvers/riccati/ilqgames/ilqgames_methods.jl")
include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("../../src/sampler/monte_carlo_sampler.jl")
include("../../src/sampler/monte_carlo_methods.jl")

include("../../src/scenarios/scenario.jl")
include("../../src/scenarios/examples/merging.jl")
include("../../src/scenarios/examples/straight.jl")
include("../../src/scenarios/examples/t_intersection.jl")

include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

include("../../src/scenarios/scenario_visualization.jl")
include("../../src/scenarios/adaptive_plot.jl")

include("../../src/utils/constraints.jl")
include("../../src/utils/monte_carlo_visualization_latex.jl")
include("../../src/utils/monte_carlo_visualization.jl")
include("../../src/utils/plot_visualization.jl")
include("../../src/utils/tests.jl")
include("../../src/utils/timing.jl")



# using ALGAMES
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization


# Define the dynamics model of the game.
struct InertialUnicycleGame{T} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
    mp::T
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px::Vector{Vector{Int}} # Indices of the each player's x and y positions
    p::Int  # Number of players
end
InertialUnicycleGame() = InertialUnicycleGame(
	12,
	6,
	1.0,
	[[1,2],[3,4],[5,6]],
	[[1,2],[5,6],[9,10]],
	3)
Base.size(::InertialUnicycleGame) = 12,6,[[1,2],[3,4],[5,6]],3 # n,m,pu,p

# Instantiate dynamics model
model = InertialUnicycleGame()
n,m,pu,p = size(model)
T = Float64
px = model.px
function TO.dynamics(model::InertialUnicycleGame, x, u) # Non memory allocating dynamics
    qd1 = @SVector [cos(x[3]), sin(x[3])]
    qd1 *= x[4]
    qd2 = @SVector [cos(x[7]), sin(x[7])]
    qd2 *= x[8]
    qd3 = @SVector [cos(x[11]), sin(x[11])]
    qd3 *= x[12]
    qdd1 = u[ @SVector [1,2] ]
    qdd2 = u[ @SVector [3,4] ]
    qdd3 = u[ @SVector [5,6] ]
    return [qd1; qdd1; qd2; qdd2; qd3; qdd3]
end

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.80, -0.05,  0.00, 0.60, #lane1
               -1.00, -0.05,  0.00, 0.60, #lane2
               -0.90, -0.30, pi/12, 0.63, #lane4
                ]
xf = @SVector [
                1.10, -0.05,  0.00, 0.60, #lane1
                0.70, -0.05,  0.00, 0.60, #lane2
                0.90, -0.05,  0.00, 0.60, #lane4
               ]

# Define a quadratic cost
diag_Q1 = @SVector [ # Player 1 state cost
    0., 1., 1., 1.,
    0., 0., 0., 0.,
    0., 0., 0., 0.]
diag_Q2 = @SVector [ # Player 2 state cost
    0., 0., 0., 0.,
    0., 1., 1., 1.,
    0., 0., 0., 0.]
diag_Q3 = @SVector [ # Player 3 state cost
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q1), # Players state costs
     0.1*Diagonal(diag_Q2),
     0.1*Diagonal(diag_Q3)]
Qf = [1.0*Diagonal(diag_Q1),
      1.0*Diagonal(diag_Q2),
      1.0*Diagonal(diag_Q3)]

# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[1]))),
     0.1*Diagonal(@SVector ones(length(pu[2]))),
     0.1*Diagonal(@SVector ones(length(pu[3]))),
    ]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
xs = SVector{n}(zeros(n))
us = SVector{m}(zeros(m))
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
car_radius = 0.08
car_radii = [car_radius for i=1:p]
actors_types = [:car, :car, :car]
road_length = 6.0
road_width = 0.30
ramp_length = 3.2
ramp_angle = pi/12
scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, car_radii, actors_types)


# Create constraint sets
algames_conSet = ConstraintSet(n,m,N)
ilqgames_conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, car_radii, px, p, con_inds)
add_collision_avoidance(ilqgames_conSet, car_radii, px, p, con_inds)

# Add scenario specific constraints (road boundaries)
add_scenario_constraints(algames_conSet, scenario, px, con_inds; constraint_type=:constraint)
add_scenario_constraints(ilqgames_conSet, scenario, px, con_inds; constraint_type=:constraint);


algames_prob = GameProblem(model, obj, algames_conSet, x0, xf, Z, N, tf);
ilqgames_prob = GameProblem(model, obj, ilqgames_conSet, x0, xf, Z, N, tf);

# AlGAMES
algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10)
algames_solver = DirectGamesSolver(algames_prob, algames_opts);

# iLQGames
ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    iterations_linesearch=10)
ilqgames_solver = PenaltyiLQGamesSolver(ilqgames_prob, ilqgames_opts)

# Tune this penalty parameter to get better constraint satisfaction for iLQGames
pen = 100.0 * ones(length(ilqgames_solver.constraints))
set_penalty!(ilqgames_solver, pen);



@time solve!(algames_solver);
@time solve!(ilqgames_solver);

reset!(algames_solver, reset_type=:full)
algames_solver.opts.log_level = TO.Logging.Warn
@btime timing_solve(algames_solver);

reset!(ilqgames_solver, reset_type=:full)
ilqgames_solver.opts.log_level = TO.Logging.Warn
@btime timing_solve(ilqgames_solver);

visualize_optimality_merit(algames_solver)
visualize_α(algames_solver)
visualize_cmax(algames_solver)

visualize_α(ilqgames_solver)
visualize_cmax(ilqgames_solver)
