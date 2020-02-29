using TrajectoryOptimization#v1.3
using StaticArrays
using Statistics
using LinearAlgebra
using BenchmarkTools
using SparseArrays
using PartedArrays
# using Plots
using Logging
const TO = TrajectoryOptimization
import TrajectoryOptimization:
    dynamics,
    RK3,
    AbstractModel,
    KnotPoint,
    Traj,
    BoundConstraint,
    GoalConstraint,
    ConstraintVals,
    ConstraintSets,
    Problem,
    ALSolver
using Parameters


include("../../src/solvers/game_model.jl")
include("../../src/solvers/game_problem.jl")

include("../../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("../../src/solvers/riccati/ilqgames/ilqgames_methods.jl")

include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("../../src/solvers/riccati/algames/algames_solver.jl")
include("../../src/solvers/riccati/algames/algames_methods.jl")

include("../../src/solvers/direct/direct_solver.jl")
include("../../src/solvers/direct/direct_methods.jl")
include("../../src/solvers/direct/direct_helpers.jl")
include("../../src/solvers/direct/direct_core.jl")

include("../../src/solvers/inds_helpers.jl")
include("../../src/solvers/cost_helpers.jl")

include("../../src/utils/constraints.jl")
include("../../src/utils/timing.jl")

include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

include("../../src/scenarios/scenario.jl")
include("../../src/scenarios/scenario_visualization.jl")
include("../../src/scenarios/old/straight_constraints.jl")
include("../../src/scenarios/old/intersection_constraints.jl")
include("../../src/scenarios/examples/merging.jl")
include("../../src/scenarios/examples/straight.jl")
include("../../src/scenarios/examples/t_intersection.jl")


# Define the dynamics model of the game.
struct DoubleIntegratorGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Vector{Vector{Int}}
    p::Int
end
DoubleIntegratorGame() = DoubleIntegratorGame(
    12, 6, 1.0, [[1,2],[3,4],[5,6]], 3)
Base.size(::DoubleIntegratorGame) = 12,6,[[1,2],[3,4],[5,6]],3

# Instantiate dynamics model
model = DoubleIntegratorGame()
n,m,pl,p = size(model)
T = Float64
plx = [[1,2], [5,6], [9,10]]

function dynamics(model::DoubleIntegratorGame, x, u)
    mp = model.mp  # mass of the point mass in kg (10)
    p = model.p  # number of players
    pl = model.pl  # control vector partition for each player
    q1 = x[ @SVector [1,2] ]
    qd1 = x[ @SVector [3,4] ]
    q2 = x[ @SVector [5,6] ]
    qd2 = x[ @SVector [7,8] ]
    q3 = x[ @SVector [9,10] ]
    qd3 = x[ @SVector [11,12] ]
    control1 = @SVector [u[pl_ind] for pl_ind in pl[1]]
    control2 = @SVector [u[pl_ind] for pl_ind in pl[2]]
    control3 = @SVector [u[pl_ind] for pl_ind in pl[3]]
    qdd1 = control1/mp
    qdd2 = control2/mp
    qdd3 = control3/mp
    return [qd1; qdd1; qd2; qdd2; qd3; qdd3]
end

# Discretization info
tf = 2.0  # final time
N = 41   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.80, -0.05,  0.90,  0.00, #lane1 player 2 goes first
               -1.00, -0.05,  0.85,  0.00, #lane2 player 2 goes first
               -0.90, -0.35,  0.85,  0.20, #lane2 player 2 goes first
                ]
xf = @SVector [
                1.10, -0.05,  0.90,  0.0, #lane1
                0.70, -0.05,  0.90,  0.0, #lane1
                0.90, -0.05,  0.90,  0.0, #lane1
               ]
dxf = @SVector [
                1.00,  0.00,  0.00,  0.00, #lane1
                1.00,  0.00,  0.00,  0.00, #lane2
                1.00,  0.00,  0.00,  0.00, #lane2
               ]



# Define a quadratic cost
diag_Q = @SVector [
    0., 10., 1., 1.,
    0., 10., 1., 1.,
    0., 10., 1., 1.]
Q = [0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q)]
Qf = [1.0*Diagonal(diag_Q),
      1.0*Diagonal(diag_Q),
      1.0*Diagonal(diag_Q)]

# Q = [0.1*Diagonal(@SVector ones(n)),
#      0.1*Diagonal(@SVector ones(n)),
#      0.1*Diagonal(@SVector ones(n))]
# Qf = [1.0*Diagonal(@SVector ones(n)),
#      1.0*Diagonal(@SVector ones(n)),
#      1.0*Diagonal(@SVector ones(n))]
R = [0.1*Diagonal(@SVector ones(length(pl[1]))),
    0.1*Diagonal(@SVector ones(length(pl[2]))),
    0.1*Diagonal(@SVector ones(length(pl[3]))),
    ]
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector zeros(m)
U0 = [TO.copy(u0) for k = 1:N]
xs = NaN*@SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
car_radius = 0.80/10.
car_radii = [car_radius for i=1:p]
road_length = 2.00
road_width = 0.30
ramp_length = 1.2
ramp_angle = pi/12
# obs_radius = 2.0/10. - car_radius
# obs = [0.0, -0.20]
scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, car_radii)

# Collision Avoidance
conSet = ConstraintSet(n,m,N)
add_collision_avoidance(conSet, car_radii, plx, p)
add_scenario_constraints(conSet, scenario, plx, N, n; constraint_type=:constraint)
# add_scenario_constraints(conSet, scenario, plx, N, n; constraint_type=:constraint)


include("../../src/solvers/direct/direct_solver.jl")
prob = GameProblem(model, obj, conSet, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=5,
    info_pattern=:open_loop)
# opts_ilqgames = iLQGamesSolverOptions{T}(
#     iterations=10,
#     info_pattern=:feedback)
# opts_algames = ALGamesSolverOptions{T}()
solver_directgames = DirectGamesSolver(prob, opts_directgames)
# solver_ilqgames = iLQGamesSolver(prob, opts_ilqgames)
# solver_algames = ALGamesSolver(prob, opts_algames)
opts_penalty = PenaltyiLQGamesSolverOptions{T}(iterations=100,
    line_search_lower_bound=-1.0,
    line_search_upper_bound=10.0)
solver_penalty = PenaltyiLQGamesSolver(prob, opts_penalty)
pen = ones(length(solver_penalty.constraints))*10000.0
set_penalty!(solver_penalty, pen)
# initial_controls!(prob,U0)
# @time solve!(solver_algames)
# @time solve!(solver_ilqgames)
@time solve!(solver_penalty)
# @time solve!(solver_directgames)


# X = TO.states(solver_directgames)
# U = TO.controls(solver_directgames)
# include("../../src/utils/plot_visualization.jl")
# using Plots
# visualize_trajectory_car(solver_directgames,pl,plx,car_radius,l1)
#
@time solve!(solver_directgames)
full_reset!(solver_directgames, prob)
@time solve!(solver_directgames)
full_reset!(solver_directgames, prob)
@profiler solve!(solver_directgames)
full_reset!(solver_directgames, prob)
@btime timing_solve(solver_directgames, prob)

# @time solve!(solver_ilqgames)
# reset!(solver_ilqgames)
# @time solve!(solver_ilqgames)
# reset!(solver_ilqgames)
# @profiler solve!(solver_ilqgames)
# reset!(solver_ilqgames)
# @btime timing_solve(solver_ilqgames)

@time solve!(solver_penalty)
reset!(solver_penalty)
@time solve!(solver_penalty)
reset!(solver_penalty)
# @profiler solve!(solver_penalty)
# reset!(solver_penalty)
@btime timing_solve(solver_penalty)


include("../../src/utils/plot_visualization.jl")
using Plots
X = TO.states(solver_penalty)
U = TO.controls(solver_penalty)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_penalty,pl,plx,car_radius,road_length)
visualize_collision_avoidance(solver_penalty)
visualize_circle_collision(solver_penalty)
visualize_boundary_collision(solver_penalty)
# visualize_dynamics(solver_penalty)
# visualize_optimality_merit(solver_penalty)
# visualize_H_cond(solver_penalty)
visualize_α(solver_penalty)
visualize_cmax(solver_penalty)


# include("../../src/utils/plot_visualization.jl")
# using Plots
# X = TO.states(solver_ilqgames)
# U = TO.controls(solver_ilqgames)
# visualize_state(X)
# visualize_control(U,pl)
# visualize_trajectory_car(solver_ilqgames,pl,plx,car_radius,l1)
# visualize_collision_avoidance(solver_ilqgames)
# visualize_circle_collision(solver_ilqgames)
# visualize_boundary_collision(solver_ilqgames)
# visualize_dynamics(solver_ilqgames)
# visualize_optimality_merit(solver_ilqgames)
# visualize_H_cond(solver_ilqgames)
# visualize_α(solver_ilqgames)
# visualize_cmax(solver_ilqgames)


X = TO.states(solver_directgames)
U = TO.controls(solver_directgames)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_directgames,pl,plx,car_radius,road_length)
visualize_collision_avoidance(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_boundary_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)

# include("../../src/utils/meshcat_visualization.jl")
# animation(solver_ilqgames, plx, car_radius, l1, l2)

include("../../src/utils/meshcat_visualization.jl")
animation(solver_penalty, scenario, plx)
animation(solver_directgames, scenario, plx)




a = 100
a = 100
a = 100
a = 100
a = 100




include("../../src/solvers/direct/direct_solver.jl")
include("../../src/solvers/direct/direct_methods.jl")

include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

solver_directgames = DirectGamesSolver(prob, opts_directgames)
initial_controls!(prob,U0)
@time solve!(solver_directgames)
full_reset!(solver_directgames, prob)
opts_mpc_gamessolver = MPCGamesSolverOptions{T}(iterations=200)
mpc_solver = MPCGamesSolver(solver_directgames, dxf, opts_mpc_gamessolver)
reset!(mpc_solver)
solve!(mpc_solver)

a = 100
a = 100
a = 100
a = 100
a = 100
