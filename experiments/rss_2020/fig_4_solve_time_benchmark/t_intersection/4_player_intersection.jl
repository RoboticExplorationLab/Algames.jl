using TrajectoryOptimization#v1.3
using StaticArrays
using Statistics
using LinearAlgebra
using BenchmarkTools
using SparseArrays
using PartedArrays
using Dates
using Random
using JLD2
using Plots
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
include("../../src/utils/plot_visualization.jl")

include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

include("../../src/scenarios/scenario.jl")
include("../../src/scenarios/scenario_visualization.jl")
include("../../src/scenarios/old/straight_constraints.jl")
include("../../src/scenarios/old/intersection_constraints.jl")
include("../../src/scenarios/examples/merging.jl")
include("../../src/scenarios/examples/straight.jl")
include("../../src/scenarios/examples/t_intersection.jl")


include("../../src/sampler/monte_carlo_sampler.jl")
include("../../src/sampler/monte_carlo_methods.jl")


# Define the dynamics model of the game.
struct InertialUnicycleGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Vector{Vector{Int}}
    p::Int
end
InertialUnicycleGame() = InertialUnicycleGame(
    16, 8, 1.0, [[1,2],[3,4],[5,6],[7,8]], 4)
Base.size(::InertialUnicycleGame) = 16,8,[[1,2],[3,4],[5,6],[7,8]],4

# Instantiate dynamics model
model = InertialUnicycleGame()
n,m,pl,p = size(model)
T = Float64
plx = [[1,2], [5,6], [9,10], [13,14]]


function dynamics(model::InertialUnicycleGame, x, u)
    qd1 = @SVector [cos(x[3]), sin(x[3])]
    qd1 *= x[4]
    qd2 = @SVector [cos(x[7]), sin(x[7])]
    qd2 *= x[8]
    qd3 = @SVector [cos(x[11]), sin(x[11])]
    qd3 *= x[12]
    qd4 = @SVector [cos(x[15]), sin(x[15])]
    qd4 *= x[16]
    qdd1 = u[ @SVector [1,2] ]
    qdd2 = u[ @SVector [3,4] ]
    qdd3 = u[ @SVector [5,6] ]
    qdd4 = u[ @SVector [7,8] ]
    return [qd1; qdd1; qd2; qdd2; qd3; qdd3; qd4; qdd4]
end


# Discretization info
tf = 3.0  # final time
N = 41   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.50, -0.15,  0.00, 0.60, #lane1
                1.40,  0.15,  pi, 0.30, #lane2
               # -1.00, -0.05,  0.00, 0.60, #lane2
               # 0.15, -0.60,  pi/2, 0.60, #lane4
               0.15, -0.30,  pi/2, 0.80, #lane4
               0.43, -0.30,  pi/2, 0.10, #lane4
                ]
xf = @SVector [
                1.30, -0.15,  0.00, 0.60, #lane1
               -0.30,  0.15,  pi, 0.30, #lane2
                # 0.70, -0.05,  0.00, 0.60, #lane2
               # -0.90,  0.15,  pi, 0.80, #lane4
              -0.60,  0.15,  pi, 0.80, #lane4
               0.43,  0.35,  pi/2, 0.30, #lane4
               ]


dxf = @SVector [
               0.60,  0.00,  0.00,  0.00, #lane1
              -0.30,  0.00,  0.00,  0.00, #lane2
               # 1.00,  0.00,  0.00,  0.00, #lane2
              -0.60,  0.00,  0.00,  0.00, #lane2
               0.00,  0.10,  0.00,  0.00, #lane2
              ]



# Define a quadratic cost
diag_Q = @SVector [
    0., 1., 1., 1.,
    0., 1., 1., 1.,
    0., 1., 1., 1.,
    0., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q)]
Qf = [1.0*Diagonal(diag_Q),
      1.0*Diagonal(diag_Q),
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
    0.1*Diagonal(@SVector ones(length(pl[2]))),
    0.1*Diagonal(@SVector ones(length(pl[3]))),
    ]
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector zeros(m)
U0 = [TO.copy(u0) for k = 1:N]
xs = @SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)


# Build problem
car_radius = 0.80/10.
car_radii = [car_radius for i=1:p]
actors_types = [:car, :car, :car, :pedestrian]
# road_length = 6.0
# road_width = 0.30
# ramp_length = 3.2
# ramp_angle = pi/12
top_road_length = 4.0
road_width = 0.60
bottom_road_length = 1.0
cross_width = 0.25
bound_radius = 0.05
# obs_radius = 2.0/10. - car_radius
# obs = [0.0, -0.20]
# lanes = [1, 2, 3, 5]
lanes = [1, 2, 3, 5]
scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length, cross_width, car_radii, actors_types, bound_radius)

# Collision Avoidance
conSet_direct = ConstraintSet(n,m,N)
conSet_penalty = ConstraintSet(n,m,N)
add_collision_avoidance(conSet_direct, car_radii, plx, p)
add_collision_avoidance(conSet_penalty, car_radii, plx, p)
add_scenario_constraints(conSet_direct, scenario, lanes, plx, N, n; constraint_type=:constraint)
add_scenario_constraints(conSet_penalty, scenario, lanes, plx, N, n; constraint_type=:constraint)


# include("../../src/solvers/direct/direct_solver.jl")
prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
prob_penalty = GameProblem(model, obj, conSet_penalty, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    info_pattern=:open_loop)
    solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
opts_penalty = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    iterations_linesearch=5,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.02)
solver_penalty = PenaltyiLQGamesSolver(prob_penalty, opts_penalty)
pen = ones(length(solver_penalty.constraints))*1000.0
set_penalty!(solver_penalty, pen)
@time solve!(solver_penalty)
# full_reset!(solver_directgames)
@time solve!(solver_directgames)
# rollout!(solver_directgames)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 200


#
@time solve!(solver_directgames)
full_reset!(solver_directgames)
@time solve!(solver_directgames)
# full_reset!(solver_directgames)
# @profiler solve!(solver_directgames)
full_reset!(solver_directgames)
direct_bench = @benchmark timing_solve(solver_directgames)
direct_mean = mean(direct_bench.times)/1e6
direct_std = sqrt(var(direct_bench.times))/1e6
length(direct_bench.times)

@time solve!(solver_penalty)
full_reset!(solver_penalty)
@time solve!(solver_penalty)
# full_reset!(solver_penalty)
# @profiler solve!(solver_penalty)
full_reset!(solver_penalty)
penalty_bench = @benchmark timing_solve(solver_penalty)
penalty_mean = mean(penalty_bench.times)/1e6
penalty_std = sqrt(var(penalty_bench.times))/1e6
length(penalty_bench.times)


include("../../src/utils/plot_visualization.jl")
using Plots

X = TO.states(solver_penalty)
U = TO.controls(solver_penalty)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_penalty,pl,plx)
visualize_collision_avoidance(solver_penalty)
visualize_circle_collision(solver_penalty)
visualize_boundary_collision(solver_penalty)
# visualize_dynamics(solver_penalty)
# visualize_optimality_merit(solver_penalty)
# visualize_H_cond(solver_penalty)
visualize_α(solver_penalty)
visualize_cmax(solver_penalty)

X = TO.states(solver_directgames)
U = TO.controls(solver_directgames)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_directgames,pl,plx)
visualize_collision_avoidance(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_boundary_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)

animation(solver_directgames, scenario, plx)
# animation(solver_penalty, scenario, plx)



a = 10
a = 10
a = 10
a = 10
a = 10


#
# include("../../src/solvers/MPC/mpc_solver.jl")
# include("../../src/solvers/MPC/mpc_methods.jl")
#
# solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
# initial_controls!(prob_direct,U0)
# @time solve!(solver_directgames)
# full_reset!(solver_directgames)
#
# opts_mpc_gamessolver = MPCGamesSolverOptions{T}(iterations=1000)
# mpc_solver = MPCGamesSolver(solver_directgames, dxf, opts_mpc_gamessolver)
#
# samples = 100
# times = zeros(0)
# cmax = zeros(samples)
# for k = 1:samples
#     reset!(mpc_solver)
#     state_noise = @SVector [
#         0.008, 0.008, 2*pi/72, 0.03, #+-10cm, +-10cm, +-5deg, +-2.5%
#         0.008, 0.008, 2*pi/72, 0.03,
#         0.008, 0.008, 2*pi/72, 0.03]
#     state_noise *= 5.0
#     selfish_inds = zeros(Int,0)
#     selfish_dx = SVector{0}([])
#     solve!(mpc_solver, state_noise, selfish_inds, selfish_dx; min_δt=0.00, display=false)
#
#     i = mpc_solver.stats.iterations
#     cmax[k] = maximum(mpc_solver.stats.cmax[1:i])
#     push!(times, mpc_solver.stats.solve_time[1:i]...)
# end
#
# times
# mean(times)
# var(times)
# sqrt(var(times))
# cmax
#
#
#
#
# a = 10
# a = 10
# a = 10
# a = 10
# a = 10

















using TrajectoryOptimization#v1.3
using StaticArrays
using Statistics
using LinearAlgebra
using BenchmarkTools
using SparseArrays
using PartedArrays
using Dates
using Random
using JLD2
using Plots
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

include("../../../src/solvers/game_model.jl")
include("../../../src/solvers/game_problem.jl")

include("../../../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("../../../src/solvers/riccati/ilqgames/ilqgames_methods.jl")

include("../../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("../../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("../../../src/solvers/riccati/algames/algames_solver.jl")
include("../../../src/solvers/riccati/algames/algames_methods.jl")

include("../../../src/solvers/direct/direct_solver.jl")
include("../../../src/solvers/direct/direct_methods.jl")
include("../../../src/solvers/direct/direct_helpers.jl")
include("../../../src/solvers/direct/direct_core.jl")

include("../../../src/solvers/inds_helpers.jl")
include("../../../src/solvers/cost_helpers.jl")

include("../../../src/utils/constraints.jl")
include("../../../src/utils/timing.jl")
include("../../../src/utils/plot_visualization.jl")

include("../../../src/solvers/MPC/mpc_solver.jl")
include("../../../src/solvers/MPC/mpc_methods.jl")

include("../../../src/scenarios/scenario.jl")
include("../../../src/scenarios/scenario_visualization.jl")
include("../../../src/scenarios/old/straight_constraints.jl")
include("../../../src/scenarios/old/intersection_constraints.jl")
include("../../../src/scenarios/examples/merging.jl")
include("../../../src/scenarios/examples/straight.jl")
include("../../../src/scenarios/examples/t_intersection.jl")


include("../../../src/sampler/monte_carlo_sampler.jl")
include("../../../src/sampler/monte_carlo_methods.jl")


# Define the dynamics model of the game.
struct InertialUnicycleGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Vector{Vector{Int}}
    p::Int
end
InertialUnicycleGame() = InertialUnicycleGame(
    16, 8, 1.0, [[1,2],[3,4],[5,6],[7,8]], 4)
Base.size(::InertialUnicycleGame) = 16,8,[[1,2],[3,4],[5,6],[7,8]],4

# Instantiate dynamics model
model = InertialUnicycleGame()
n,m,pl,p = size(model)
T = Float64
plx = [[1,2], [5,6], [9,10], [13,14]]


function dynamics(model::InertialUnicycleGame, x, u)
    qd1 = @SVector [cos(x[3]), sin(x[3])]
    qd1 *= x[4]
    qd2 = @SVector [cos(x[7]), sin(x[7])]
    qd2 *= x[8]
    qd3 = @SVector [cos(x[11]), sin(x[11])]
    qd3 *= x[12]
    qd4 = @SVector [cos(x[15]), sin(x[15])]
    qd4 *= x[16]
    qdd1 = u[ @SVector [1,2] ]
    qdd2 = u[ @SVector [3,4] ]
    qdd3 = u[ @SVector [5,6] ]
    qdd4 = u[ @SVector [7,8] ]
    return [qd1; qdd1; qd2; qdd2; qd3; qdd3; qd4; qdd4]
end



# Discretization info
tf = 3.0  # final time
N = 41   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.80, -0.10,  0.00, 0.62, #lane1
               -1.00, -0.10,  0.00, 0.60, #lane2
               -0.90, -0.35, pi/12, 0.63, #lane4
               -0.90,  0.05,  0.00, 0.79, #lane4
                ]
xf = @SVector [
                1.10, -0.10,  0.00, 0.62, #lane1
                0.70, -0.10,  0.00, 0.60, #lane2
                0.90, -0.10,  0.00, 0.60, #lane4
                1.30, -0.10,  0.00, 0.79, #lane4
               ]


dxf = @SVector [
               1.00,  0.00,  0.00,  0.00, #lane1
               1.00,  0.00,  0.00,  0.00, #lane2
               1.00,  0.00,  0.00,  0.00, #lane2
               1.00,  0.00,  0.00,  0.00, #lane2
              ]



# Define a quadratic cost
diag_Q = @SVector [
    0., 10., 1., 1.,
    0., 10., 1., 1.,
    0., 10., 1., 1.,
    0., 10., 1., 1.]
Q = [0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q)]
Qf = [1.0*Diagonal(diag_Q),
      1.0*Diagonal(diag_Q),
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
    0.1*Diagonal(@SVector ones(length(pl[2]))),
    0.1*Diagonal(@SVector ones(length(pl[3]))),
    ]
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector zeros(m)
U0 = [TO.copy(u0) for k = 1:N]
xs = @SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)


# Build problem
car_radius = 0.80/10.
car_radii = [car_radius for i=1:p]
actors_types = [:car, :car, :car, :car]
road_length = 6.0
road_width = 0.40
ramp_length = 3.2
ramp_angle = pi/12
# obs_radius = 2.0/10. - car_radius
# obs = [0.0, -0.20]
scenario = MergingScenario(road_length, road_width, ramp_length, ramp_angle, car_radii, actors_types)

# Collision Avoidance
conSet_direct = ConstraintSet(n,m,N)
conSet_penalty = ConstraintSet(n,m,N)
add_collision_avoidance(conSet_direct, car_radii, plx, p)
add_collision_avoidance(conSet_penalty, car_radii, plx, p)
add_scenario_constraints(conSet_direct, scenario, plx, N, n; constraint_type=:constraint)
add_scenario_constraints(conSet_penalty, scenario, plx, N, n; constraint_type=:constraint)


# include("../../../src/solvers/direct/direct_solver.jl")
prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
prob_penalty = GameProblem(model, obj, conSet_penalty, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    info_pattern=:open_loop)
    solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
opts_penalty = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    iterations_linesearch=5,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.02)
solver_penalty = PenaltyiLQGamesSolver(prob_penalty, opts_penalty)
pen = ones(length(solver_penalty.constraints))*100.0
set_penalty!(solver_penalty, pen)
full_reset!(solver_penalty)
@time solve!(solver_penalty)
full_reset!(solver_directgames)
@time solve!(solver_directgames)
# rollout!(solver_directgames)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 200



@time solve!(solver_directgames)
full_reset!(solver_directgames)
@time solve!(solver_directgames)
# full_reset!(solver_directgames)
# @profiler solve!(solver_directgames)
full_reset!(solver_directgames)
direct_bench = @benchmark timing_solve(solver_directgames)
direct_mean = mean(direct_bench.times)/1e6
direct_std = sqrt(var(direct_bench.times))/1e6
length(direct_bench.times)

@time solve!(solver_penalty)
full_reset!(solver_penalty)
@time solve!(solver_penalty)
# full_reset!(solver_penalty)
# @profiler solve!(solver_penalty)
full_reset!(solver_penalty)
penalty_bench = @benchmark timing_solve(solver_penalty)
penalty_mean = mean(penalty_bench.times)/1e6
penalty_std = sqrt(var(penalty_bench.times))/1e6
length(penalty_bench.times)




include("../../../src/utils/plot_visualization.jl")
using Plots

X = TO.states(solver_penalty)
U = TO.controls(solver_penalty)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_penalty,pl,plx)
visualize_collision_avoidance(solver_penalty)
visualize_circle_collision(solver_penalty)
visualize_boundary_collision(solver_penalty)
# visualize_dynamics(solver_penalty)
# visualize_optimality_merit(solver_penalty)
# visualize_H_cond(solver_penalty)
visualize_α(solver_penalty)
visualize_cmax(solver_penalty)

X = TO.states(solver_directgames)
U = TO.controls(solver_directgames)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_directgames,pl,plx)
visualize_collision_avoidance(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_boundary_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)

animation(solver_directgames, scenario, plx)
animation(solver_penalty, scenario, plx)

a = 100
a = 100
a = 100
a = 100
a = 100

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
@btime solver_directgames.H_ \ solver_directgames.g_




# include("../../../src/solvers/MPC/mpc_solver.jl")
# include("../../../src/solvers/MPC/mpc_methods.jl")
#
# solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
# initial_controls!(prob_direct,U0)
# @time solve!(solver_directgames)
# full_reset!(solver_directgames)
# opts_mpc_gamessolver = MPCGamesSolverOptions{T}(iterations=3000)
#
#
# samples = 100
# times = zeros(0)
# cmax = zeros(samples)
# for k = 1:samples
#     reset!(mpc_solver)
#     state_noise = @SVector [
#         0.008, 0.008, 2*pi/72, 0.03, #+-10cm, +-10cm, +-5deg, +-2.5%
#         0.008, 0.008, 2*pi/72, 0.03,
#         0.008, 0.008, 2*pi/72, 0.03]
#     state_noise *= 5.0
#     selfish_inds = zeros(Int,0)
#     selfish_dx = SVector{0}([])
#     solve!(mpc_solver, state_noise, selfish_inds, selfish_dx; min_δt=0.00, display=false)
#
#     i = mpc_solver.stats.iterations
#     cmax[k] = maximum(mpc_solver.stats.cmax[1:i])
#     push!(times, mpc_solver.stats.solve_time[1:i]...)
# end
#
# times
# mean(times)
# var(times)
# sqrt(var(times))
# cmax
#








#
#
# mpc_solver = MPCGamesSolver(solver_directgames, dxf, opts_mpc_gamessolver)
# reset!(mpc_solver)
# state_noise = @SVector [
#     0.008, 0.008, 2*pi/72, 0.03, #+-10cm, +-10cm, +-5deg, +-2.5%
#     0.008, 0.008, 2*pi/72, 0.03,
#     0.008, 0.008, 2*pi/72, 0.03]
# state_noise *= 5.0
# solve!(mpc_solver, state_noise; min_δt=1e-23)
#
# i = mpc_solver.stats.iterations
# mean(mpc_solver.stats.solve_time[1:i])
# maximum(mpc_solver.stats.solve_time[1:i])
# sqrt(var(mpc_solver.stats.solve_time[1:i]))
# maximum(mpc_solver.stats.cmax[1:i])
#
# 1/ 0.00805

a = 10
a = 10
a = 10
a = 10
a = 10


















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

include("../../../../src/solvers/game_model.jl")
include("../../../../src/solvers/game_problem.jl")
include("../../../../src/solvers/cost_helpers.jl")

include("../../../../src/solvers/direct/direct_helpers.jl")

include("../../../../src/solvers/direct/direct_solver.jl")
include("../../../../src/solvers/direct/direct_methods.jl")
include("../../../../src/solvers/direct/direct_core.jl")
include("../../../../src/solvers/direct/newton_gradient.jl")
include("../../../../src/solvers/direct/newton_hessian.jl")
include("../../../../src/solvers/inds_helpers.jl")


include("../../../../src/solvers/riccati/algames/algames_solver.jl")
include("../../../../src/solvers/riccati/algames/algames_methods.jl")
include("../../../../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("../../../../src/solvers/riccati/ilqgames/ilqgames_methods.jl")
include("../../../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("../../../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("../../../../src/sampler/monte_carlo_sampler.jl")
include("../../../../src/sampler/monte_carlo_methods.jl")

include("../../../../src/scenarios/scenario.jl")
include("../../../../src/scenarios/examples/merging.jl")
include("../../../../src/scenarios/examples/straight.jl")
include("../../../../src/scenarios/examples/t_intersection.jl")

include("../../../../src/solvers/MPC/mpc_solver.jl")
include("../../../../src/solvers/MPC/mpc_methods.jl")

include("../../../../src/scenarios/scenario_visualization.jl")
include("../../../../src/scenarios/adaptive_plot.jl")

include("../../../../src/utils/constraints.jl")
include("../../../../src/utils/monte_carlo_visualization_latex.jl")
include("../../../../src/utils/monte_carlo_visualization.jl")
include("../../../../src/utils/plot_visualization.jl")
include("../../../../src/utils/tests.jl")
include("../../../../src/utils/timing.jl")



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
	16,
	8,
	1.0,
	[[1,2],[3,4],[5,6],[7,8]],
	[[1,2],[5,6],[9,10],[13,14]],
	4)
Base.size(::InertialUnicycleGame) = 16,8,[[1,2],[3,4],[5,6],[7,8]],4 # n,m,pu,p

# Instantiate dynamics model
model = InertialUnicycleGame()
n,m,pu,p = size(model)
T = Float64
px = model.px
function TO.dynamics(model::InertialUnicycleGame, x, u)
    qd1 = @SVector [cos(x[3]), sin(x[3])]
    qd1 *= x[4]
    qd2 = @SVector [cos(x[7]), sin(x[7])]
    qd2 *= x[8]
    qd3 = @SVector [cos(x[11]), sin(x[11])]
    qd3 *= x[12]
    qd4 = @SVector [cos(x[15]), sin(x[15])]
    qd4 *= x[16]
    qdd1 = u[ @SVector [1,2] ]
    qdd2 = u[ @SVector [3,4] ]
    qdd3 = u[ @SVector [5,6] ]
    qdd4 = u[ @SVector [7,8] ]
    return [qd1; qdd1; qd2; qdd2; qd3; qdd3; qd4; qdd4]
end

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.50, -0.15,  0.00, 0.60, #lane1
                1.40,  0.15,  pi, 0.30,   #lane2
                0.15, -0.30,  pi/2, 0.80, #lane4
                0.43, -0.30,  pi/2, 0.10, #lane4
                ]
xf = @SVector [
                1.30, -0.15,  0.00, 0.60, #lane1
               -0.30,  0.15,  pi, 0.30,   #lane2
               -0.60,  0.15,  pi, 0.80,   #lane4
                0.43,  0.35,  pi/2, 0.30, #lane4
               ]

# Define a quadratic cost
diag_Q1 = @SVector [ # Player 1 state cost
    0., 1., 1., 1.,
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 0., 0., 0.]
diag_Q2 = @SVector [ # Player 2 state cost
    0., 0., 0., 0.,
    0., 1., 1., 1.,
    0., 0., 0., 0.,
    0., 0., 0., 0.]
diag_Q3 = @SVector [ # Player 3 state cost
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 1., 1., 1.,
    0., 0., 0., 0.]
diag_Q4 = @SVector [ # Player 4 state cost
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q1), # Players state costs
     0.1*Diagonal(diag_Q2),
     0.1*Diagonal(diag_Q3),
     0.1*Diagonal(diag_Q4)]
Qf = [1.0*Diagonal(diag_Q1),
      1.0*Diagonal(diag_Q2),
      1.0*Diagonal(diag_Q3),
      1.0*Diagonal(diag_Q4)]

# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[1]))),
     0.1*Diagonal(@SVector ones(length(pu[2]))),
     0.1*Diagonal(@SVector ones(length(pu[3]))),
     0.1*Diagonal(@SVector ones(length(pu[4]))),
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
actors_types = [:car, :car, :car, :pedestrian]
top_road_length = 4.0
road_width = 0.60
bottom_road_length = 1.0
cross_width = 0.25
bound_radius = 0.05
lanes = [1, 2, 3, 5]
scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length, cross_width, actors_radii, actors_types, bound_radius)


# Create constraints
algames_conSet = ConstraintSet(n,m,N)
ilqgames_conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px, p, con_inds)
add_collision_avoidance(ilqgames_conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, scenario, lanes, px, con_inds; constraint_type=:constraint)
add_scenario_constraints(ilqgames_conSet, scenario, lanes, px, con_inds; constraint_type=:constraint)

algames_prob = GameProblem(model, obj, algames_conSet, x0, xf, Z, N, tf)
ilqgames_prob = GameProblem(model, obj, ilqgames_conSet, x0, xf, Z, N, tf)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=0,
    log_level=TO.Logging.Warn)
algames_solver = DirectGamesSolver(algames_prob, algames_opts)
ilqgames_opts = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.02,
    log_level=TO.Logging.Warn,
    )
ilqgames_solver = PenaltyiLQGamesSolver(ilqgames_prob, ilqgames_opts)
pen = ones(length(ilqgames_solver.constraints))*1000.0
set_penalty!(ilqgames_solver, pen);

@time timing_solve(algames_solver)
@time timing_solve(ilqgames_solver)

# @btime timing_solve(algames_solver)
# @btime timing_solve(ilqgames_solver)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 20

reset!(algames_solver, reset_type=:full)
@time solve!(algames_solver)
reset!(algames_solver, reset_type=:full)
algames_bench = @benchmark timing_solve(algames_solver)
# Mean time in ms
algames_mean = mean(algames_bench.times)/1e6
# Standard deviation in ms
algames_std = sqrt(var(algames_bench.times))/1e6

reset!(ilqgames_solver, reset_type=:full)
@time solve!(ilqgames_solver)
reset!(ilqgames_solver, reset_type=:full)
ilqgames_bench = @benchmark timing_solve(ilqgames_solver)
# Mean time in ms
ilqgames_mean = mean(ilqgames_bench.times)/1e6
# Standard deviation in ms
ilqgames_std = sqrt(var(ilqgames_bench.times))/1e6
