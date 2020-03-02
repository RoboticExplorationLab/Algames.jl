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


include("../src/solvers/game_model.jl")
include("../src/solvers/game_problem.jl")

include("../src/solvers/direct/direct_solver.jl")
include("../src/solvers/direct/direct_methods.jl")
include("../src/solvers/direct/direct_helpers.jl")
include("../src/solvers/direct/direct_core.jl")
include("../src/solvers/direct/direct_constraints.jl")
include("../src/solvers/inds_helpers.jl")
include("../src/solvers/cost_helpers.jl")
include("../src/utils/intersection_constraints.jl")
include("../src/utils/timing.jl")
include("../src/utils/scenario.jl")
include("../src/utils/straight_constraints.jl")


include("../src/solvers/MPC/mpc_solver.jl")
include("../src/solvers/MPC/mpc_methods.jl")



# Define the dynamics model of the game.
struct DoubleIntegratorGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Vector{Vector{Int}}
    p::Int
end
DoubleIntegratorGame() = DoubleIntegratorGame(
    8, 4, 1.0, [[1,2],[3,4]], 2)
Base.size(::DoubleIntegratorGame) = 8,4,[[1,2],[3,4]],2

# Instantiate dynamics model
model = DoubleIntegratorGame()
n,m,pl,p = size(model)
T = Float64
plx = [[1,2], [5,6]]

function dynamics(model::DoubleIntegratorGame, x, u)
    mp = model.mp  # mass of the point mass in kg (10)
    p = model.p  # number of players
    pl = model.pl  # control vector partition for each player
    q1 = x[ @SVector [1,2] ]
    qd1 = x[ @SVector [3,4] ]
    q2 = x[ @SVector [5,6] ]
    qd2 = x[ @SVector [7,8] ]
    control1 = @SVector [u[pl_ind] for pl_ind in pl[1]]
    control2 = @SVector [u[pl_ind] for pl_ind in pl[2]]
    qdd1 = control1/mp
    qdd2 = control2/mp
    return [qd1; qdd1; qd2; qdd2]
end

# Discretization info
tf = 2.0  # final time
N = 21   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               # -1.00,  0.10,  0.90,  0.00, #lane1 player 1 goes first
               # -1.00, -0.10,  1.10,  0.00, #lane2 player 1 goes first
               -1.00,  0.10,  1.10,  0.00, #lane1 player 2 goes first
               -1.00, -0.10,  0.90,  0.00, #lane2 player 2 goes first
               # -1.00,  0.10,  1.10,  0.00, #lane1 player 2 goes first
               #  1.00, -0.10, -0.90,  0.00, #lane2 player 2 goes first
                ]
xf = @SVector [
                1.00,  0.10,  1.10,  0.0, #lane1
                # -1.00, -0.10, -0.90,  0.0, #lane1
                1.00, -0.10,  0.90,  0.0, #lane1
               ]

dxf = @SVector [
                1.00,  0.00,  0.00,  0.00, #lane1
                1.00,  0.00,  0.00,  0.00, #lane2
               ]



# Define a quadratic cost
diag_Q = @SVector [0., 1., 1., 1., 0., 1., 1., 1.]
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
    ]
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector [0., 0., 0., 0.]
U0 = [TO.copy(u0) for k = 1:N]
xs = NaN*@SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
conSet_direct = ConstraintSet(n,m,N)
# Collision Avoidance
car_radius = 0.8/10.
car_radii = [car_radius for i=1:p]
for i = 1:p
    for j = 1:i-1
        radiusi = @SVector fill(car_radius, 1)
        radiusj = @SVector fill(car_radius, 1)
        col_con = CollisionConstraint(n, radiusi, radiusi,
            plx[i][1], plx[i][2], plx[j][1], plx[j][2])
        add_constraint!(conSet_direct, col_con, 1:N)
    end
end

# Road boundaries
l1 = 0.20
l2 = 1.0
obs_radius = 0.4/10.
obs = [0.50, -0.12]
scenario = StraightScenario(l1, l2, car_radii, obs_radius, obs)
add_scenario_constraints(conSet_direct, scenario, plx, N, n)



include("../src/solvers/direct/direct_solver.jl")
include("../src/solvers/direct/direct_methods.jl")
include("../src/solvers/direct/direct_helpers.jl")
include("../src/solvers/direct/direct_core.jl")
include("../src/solvers/direct/direct_constraints.jl")

prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{Float64}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=5,
    info_pattern=:open_loop)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
initial_controls!(prob_direct,U0)


# X = TO.states(solver_directgames)
# U = TO.controls(solver_directgames)
# include("../src/utils/plot_visualization.jl")
# using Plots
# visualize_trajectory_car(solver_directgames,pl,plx,car_radius,l1)

@time solve!(solver_directgames)
full_reset!(solver_directgames, prob_direct)
@time solve!(solver_directgames)
full_reset!(solver_directgames, prob_direct)
@profiler solve!(solver_directgames)
full_reset!(solver_directgames, prob_direct)
@btime timing_solve(solver_directgames, prob_direct)



include("../src/utils/plot_visualization.jl")
using Plots
X = TO.states(solver_directgames)
U = TO.controls(solver_directgames)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_directgames,pl,plx,car_radius,l1)
visualize_collision_avoidance(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_boundary_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)

include("../src/utils/meshcat_visualization.jl")
animation(solver_directgames, plx, car_radius, l1, l2)


include("../src/solvers/direct/direct_solver.jl")
include("../src/solvers/direct/direct_methods.jl")


include("../src/solvers/MPC/mpc_solver.jl")
include("../src/solvers/MPC/mpc_methods.jl")



solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
initial_controls!(prob_direct,U0)
@time solve!(solver_directgames)
full_reset!(solver_directgames, prob_direct)
opts_mpc_gamessolver = MPCGamesSolverOptions{T}(iterations=200)
mpc_solver = MPCGamesSolver(solver_directgames, dxf, opts_mpc_gamessolver)
reset!(mpc_solver)
solve!(mpc_solver)

a = 100


# mpc_solver = MPCGamesSolver{T,i,n,m,L}(
#     mpc_solver.solver,
#     mpc_solver.opts,
#     mpc_solver.stats,
#     mpc_solver.x0,
#     mpc_solver.xf,
#     mpc_solver.tf,
#     mpc_solver.Z,
#     mpc_solver.logger)
#
#
# mpc_solver = MPCGamesSolver{T,I,n,m,n+m}(
#     mpc_solver.solver,
#     mpc_solver.opts,
#     mpc_solver.stats,
#     mpc_solver.x0,
#     mpc_solver.xf,
#     mpc_solver.tf,
#     mpc_solver.Z,
#     mpc_solver.logger)
# solver_directgames = DirectGamesSolver(solver_directgames, sx1)
# ful_reset!(solver_directgames)
# @time solve!(solver_directgames)

# function mpc_solve!(solver::DirectGamesSolver{T}, S) where T
#     n,m,N = size(solver)
#     X = []
#     for s = 1:S
#         @time solve!(solver)
#         X = TO.states(solver_directgames)
#         visualize_trajectory_car(solver,pl,plx,car_radius,l1)
#         new_x0 = SVector{n}(randn(n)*0.005 + discrete_dynamics(solver.model, solver.Z[1]))
#         solver = DirectGamesSolver(solver, new_x0)
#         mpc_reset!(solver)
#         # @show new_x0
#
#         push!(X, new_x0)
#     end
#     return X
# end

full_reset!(solver_directgames, prob_direct)
X = mpc_solve!(solver_directgames, 30)
# for k = 1:N
#     ix = solver_directgames.Z[k]._x
#     @show typeof(X[k])
#     @show typeof(solver_directgames.Z[k].z[ix])
#     solver_directgames.Z[k].z[ix] .= X[k]
# end

a = 100
a = 100
a = 100
a = 100
a = 100


full_reset!(solver_directgames, prob_direct)
@time solve!(solver_directgames)
solver_directgames.obj[1].cost[1].Q *= 20
full_reset!(solver_directgames, prob_direct)
@time solve!(solver_directgames)
