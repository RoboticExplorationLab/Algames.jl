
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
struct DoubleIntegratorGame{T} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
    mp::T   # Mass of the point mass double integrator
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px::Vector{Vector{Int}} # Indices of the each player's x and y positions
    p::Int  # Number of players
end
DoubleIntegratorGame() = DoubleIntegratorGame(
	8,
	4,
	1.0,
	[[1,2],[3,4]],
	[[1,2],[5,6]],
	2)
Base.size(::DoubleIntegratorGame) = 8,4,[[1,2],[3,4]],2 # n,m,pu,p

# Instantiate dynamics model
model = DoubleIntegratorGame()
n,m,pu,p = size(model)
T = Float64
px = model.px
function TO.dynamics(model::DoubleIntegratorGame, x, u) # Non memory allocating dynamics
	mp = model.mp  # mass of the point mass in kg (10)
    p = model.p  # number of players
    pu = model.pu  # control vector partition for each player
    q1 = x[ @SVector [1,2] ]
    qd1 = x[ @SVector [3,4] ]
    q2 = x[ @SVector [5,6] ]
    qd2 = x[ @SVector [7,8] ]
    control1 = @SVector [u[pu_ind] for pu_ind in pu[1]]
    control2 = @SVector [u[pu_ind] for pu_ind in pu[2]]
    qdd1 = control1/mp
    qdd2 = control2/mp
    return [qd1; qdd1; qd2; qdd2]
end

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
			   -0.50,  0.10,  0.50,  0.00, #player 1
			   -0.50, -0.10,  0.40,  0.00, #player 2
                ]
xf = @SVector [
                0.50, -0.10,  0.40,  0.00, # player 1
                0.50,  0.10,  0.50,  0.80, # player 2
               ]

# Define a quadratic cost
diag_Q1 = @SVector [ # Player 1 state cost
    1., 1., 1., 1.,
    0., 0., 0., 0.]
diag_Q2 = @SVector [ # Player 2 state cost
    0., 0., 0., 0.,
    1., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q1), # Players state costs
     0.1*Diagonal(diag_Q2)]
Qf = [1.0*Diagonal(diag_Q1),
      1.0*Diagonal(diag_Q2)]

# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[1]))),
     0.1*Diagonal(@SVector ones(length(pu[2]))),
    ]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
xs = SVector{n}(zeros(n))
us = SVector{m}(zeros(m))
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
actors_types = [:car, :car]

# Create constraints
conSet_direct = ConstraintSet(n,m,N)
conSet_penalty = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(conSet_direct, actors_radii, px, p, con_inds)
add_collision_avoidance(conSet_penalty, actors_radii, px, p, con_inds)

u_min = - SVector{m}(ones(m))
u_max = + SVector{m}(ones(m))
con = BoundConstraint(n,m,u_min=u_min,u_max=u_max)
add_constraint!(conSet_direct, con, con_inds)
prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    info_pattern=:open_loop,
	optimality_constraint_tolerance=1e-2,
	log_level=Logging.Debug)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
reset!(solver_directgames, reset_type=:full)
@time solve!(solver_directgames)

@btime timing_solve(solver_directgames)

# reset!(solver_directgames, reset_type=:full)
# @profiler GS.solve!(solver_directgames)

# X = TO.states(solver_directgames)
# U = TO.controls(solver_directgames)
# visualize_state(X)
# visualize_control(U,pu)
visualize_trajectory_car(solver_directgames)
visualize_collision_avoidance(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_boundary_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)
