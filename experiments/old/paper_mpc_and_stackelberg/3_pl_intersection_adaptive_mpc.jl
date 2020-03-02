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
include("../../src/scenarios/scenario_mpc_visualization.jl")
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
    12, 6, 1.0, [[1,2],[3,4],[5,6]], 3)
Base.size(::InertialUnicycleGame) = 12,6,[[1,2],[3,4],[5,6]],3

# Instantiate dynamics model
model = InertialUnicycleGame()
n,m,pl,p = size(model)
T = Float64
plx = [[1,2], [5,6], [9,10]]


function dynamics(model::InertialUnicycleGame, x, u)
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
N = 41   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.50, -0.15,  0.00, 0.60, #lane1
                1.40,  0.15,  pi, 0.60, #lane2
               # -1.00, -0.05,  0.00, 0.60, #lane2
               # 0.15, -0.60,  pi/2, 0.60, #lane4
               # 0.15, -0.30,  pi/2, 0.80, #lane4
               0.43, -0.30,  pi/2, 0.10, #lane4
                ]
xf = @SVector [
                1.30, -0.15,  0.00, 0.60, #lane1
               -0.30,  0.15,  pi, 0.60, #lane2
                # 0.70, -0.05,  0.00, 0.60, #lane2
               # -0.90,  0.15,  pi, 0.80, #lane4
              # -0.60,  0.15,  pi, 0.80, #lane4
               0.43,  0.35,  pi/2, 0.30, #lane4
               ]


dxf = @SVector [
               0.60,  0.00,  0.00,  0.00, #lane1
              -0.60,  0.00,  0.00,  0.00, #lane2
               # 1.00,  0.00,  0.00,  0.00, #lane2
              # -0.60,  0.00,  0.00,  0.00, #lane2
               0.00,  0.10,  0.00,  0.00, #lane2
              ]



# Define a quadratic cost
diag_Q = @SVector [
    0., 1., 1., 1.,
    0., 1., 1., 1.,
    # 0., 1., 1., 1.,
    0., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q),
     # 0.1*Diagonal(diag_Q),
     0.1*Diagonal(diag_Q)]
Qf = [1.0*Diagonal(diag_Q),
      1.0*Diagonal(diag_Q),
      # 1.0*Diagonal(diag_Q),
      1.0*Diagonal(diag_Q)]

# Q = [0.1*Diagonal(@SVector ones(n)),
#      0.1*Diagonal(@SVector ones(n)),
#      0.1*Diagonal(@SVector ones(n))]
# Qf = [1.0*Diagonal(@SVector ones(n)),
#      1.0*Diagonal(@SVector ones(n)),
#      1.0*Diagonal(@SVector ones(n))]
R = [0.1*Diagonal(@SVector ones(length(pl[1]))),
    0.1*Diagonal(@SVector ones(length(pl[2]))),
    # 0.1*Diagonal(@SVector ones(length(pl[2]))),
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
actors_types = [:car, :car, :pedestrian]
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
lanes = [1, 2, 5]
scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length,
    cross_width, car_radii, actors_types, bound_radius)

# Collision Avoidance
conSet_direct = ConstraintSet(n,m,N)
conSet_penalty = ConstraintSet(n,m,N)
add_collision_avoidance(conSet_direct, car_radii, plx, p)
add_collision_avoidance(conSet_penalty, car_radii, plx, p)
add_scenario_constraints(conSet_direct, scenario, lanes, plx, N, n; constraint_type=:constraint)
add_scenario_constraints(conSet_penalty, scenario, lanes, plx, N, n; constraint_type=:constraint)


include("../../src/solvers/direct/direct_solver.jl")
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
    gradient_norm_tolerance=1e-1,
    cost_tolerance=1e-4,
    iterations_linesearch=5,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05)
solver_penalty = PenaltyiLQGamesSolver(prob_penalty, opts_penalty)
pen = ones(length(solver_penalty.constraints))*1000.0
set_penalty!(solver_penalty, pen)
@time solve!(solver_penalty)
# full_reset!(solver_directgames)
@time solve!(solver_directgames)
# rollout!(solver_directgames)
still_animation(solver_directgames, scenario, plx)
# @time solve!(solver_directgames)
# full_reset!(solver_directgames)
# @time solve!(solver_directgames)
# full_reset!(solver_directgames)
# # @profiler solve!(solver_directgames)
# full_reset!(solver_directgames)
# @btime timing_solve(solver_directgames)


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
animation(solver_penalty, scenario, plx)


include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
initial_controls!(prob_direct,U0)
@time solve!(solver_directgames)
full_reset!(solver_directgames)
opts_mpc_gamessolver = MPCGamesSolverOptions{T}(iterations=3000, mpc_tf=5.0)
mpc_solver = MPCGamesSolver(solver_directgames, dxf, opts_mpc_gamessolver)
reset!(mpc_solver)
state_noise = @SVector [
    0.008, 0.008, 2*pi/72, 0.03, #+-10cm, +-10cm, +-5deg, +-2.5%
    0.008, 0.008, 2*pi/72, 0.03,
    0.008, 0.008, 2*pi/72, 0.03]
state_noise *= 5.0

selfish_inds = [9,10,11,12]
selfish_dx = @SVector [0., 0.1, 0., 0.]

# selfish_inds = zeros(Int,0)
# selfish_dx = SVector{0}([])
solve!(mpc_solver, state_noise, selfish_inds,
    selfish_dx; min_δt=0.02, display=true)

# animation(mpc_solver, scenario, plx)



prob_direct_clean = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
opts_directgames_clean = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    info_pattern=:open_loop)
solver_directgames_clean = DirectGamesSolver(prob_direct_clean, opts_directgames_clean)
# full_reset!(solver_directgames)
@time solve!(solver_directgames_clean)


function plot_velocity(solver_mpc, solver)
    plt = plot()
    idx = 4
    iter = solver_mpc.stats.iterations
    vel_real = [TO.states(solver_mpc)[k][idx] for k=1:iter]
    vel_pred = [TO.states(solver)[k][idx] for k=1:N]
    max_vel = max(maximum(vel_pred), maximum(vel_real))
    vel_real ./= max_vel
    vel_pred ./= max_vel


    time_real = cumsum(solver_mpc.stats.solve_time[1:iter]).- solver_mpc.stats.solve_time[1]
    time_pred = [(k-1)*solver.tf/(N-1) for k=1:N]
    max_time = max(maximum(time_pred), maximum(time_real))
    time_real ./= max_time
    time_pred ./= max_time
    plot!(time_pred, vel_pred,
        # color=:blue,
        # marker=:circle,
        linewidth=3.0,
        linestyle=:dot,
        title="Velocity",
        label="vel_pred")
    plot!(time_real, vel_real,
        # color=:blue,
        # marker=:circle,
        linewidth=3.0,
        linestyle=:dot,
        title="Velocity",
        label="vel_real")
    display(plt)
    return nothing
end

plot_velocity(mpc_solver, solver_directgames_clean)


using PGFPlotsX
function latex_plot_velocity(solver_mpc, solver)
    idx = 4
    iter = solver_mpc.stats.iterations
    vel_real = [TO.states(solver_mpc)[k][idx] for k=1:iter]
    vel_pred = [TO.states(solver)[k][idx] for k=1:N]
    max_vel = max(maximum(vel_pred), maximum(vel_real))
    vel_real ./= max_vel
    vel_pred ./= max_vel


    time_real = cumsum(solver_mpc.stats.solve_time[1:iter]).- solver_mpc.stats.solve_time[1]
    time_pred = [(k-1)*solver.tf/(N-1) for k=1:N]
    max_time = max(maximum(time_pred), maximum(time_real))
    time_real ./= max_time
    time_pred ./= max_time

    x = range(-1; stop = 1, length = 51) # so that it contains 1/0
    axis = @pgf Axis(
        {
            "legend pos=south east",
            ymajorgrids,
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
            "x label style={at={(axis description cs:0.5,-0.20)},anchor=north}",
            "y label style={at={(axis description cs:-0.10,0.5)},rotate=0,anchor=south}",
            "xlabel={Scaled Time}",
            "ylabel={Scaled Velocity}",
            xmajorgrids = false,
            xmin = 0.00,   xmax = 1.00,
            ymin = 0.00, #  ymax = 1.00,

        },
        Plot(
            {
            "thick",
            "orange",
            no_marks,
            },
            Coordinates(time_pred, vel_pred)
        ),
        LegendEntry("Initial Plan"),
        Plot(
            {
            "thick",
            "blue",
            no_marks,
            },
            Coordinates(time_real, vel_real)
        ),
        LegendEntry("MPC"),
    )
    pgfsave("plots/tikz/velocity_"*string(solver_mpc.stats.iterations)*".tikz",
        axis; include_preamble=false, dpi = 600)
    return axis
end

latex_plot_velocity(mpc_solver, solver_directgames_clean)




a = 10
a = 10
a = 10
a = 10
a = 10


# i = mpc_solver.stats.iterations
# mean(mpc_solver.stats.solve_time[1:i])
# maximum(mpc_solver.stats.solve_time[1:i])
# sqrt(var(mpc_solver.stats.solve_time[1:i]))
# maximum(mpc_solver.stats.cmax[1:i])
