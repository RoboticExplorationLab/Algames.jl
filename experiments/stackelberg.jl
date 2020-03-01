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

include("../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("../src/solvers/riccati/ilqgames/ilqgames_methods.jl")

include("../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("../src/solvers/riccati/algames/algames_solver.jl")
include("../src/solvers/riccati/algames/algames_methods.jl")

include("../src/solvers/direct/direct_solver.jl")
include("../src/solvers/direct/direct_methods.jl")
include("../src/solvers/direct/direct_helpers.jl")
include("../src/solvers/direct/direct_core.jl")

include("../src/solvers/inds_helpers.jl")
include("../src/solvers/cost_helpers.jl")

include("../src/utils/constraints.jl")
include("../src/utils/timing.jl")

include("../src/solvers/MPC/mpc_solver.jl")
include("../src/solvers/MPC/mpc_methods.jl")

include("../src/scenarios/scenario.jl")
include("../src/scenarios/scenario_visualization.jl")
include("../src/scenarios/old/straight_constraints.jl")
include("../src/scenarios/old/intersection_constraints.jl")
include("../src/scenarios/examples/merging.jl")
include("../src/scenarios/examples/straight.jl")
include("../src/scenarios/examples/t_intersection.jl")


include("../src/sampler/monte_carlo_sampler.jl")
include("../src/sampler/monte_carlo_methods.jl")


# Define the dynamics model of the game.
struct DoubleIntegrator{T} <: AbstractModel
    n::Int
    m::Int
    mp::T
    p::Int
end
DoubleIntegrator() = DoubleIntegrator(
    4, 2, 1.0, 1)
Base.size(::DoubleIntegrator) = 4, 2

# Instantiate dynamics model
model = DoubleIntegrator()
n,m = size(model)
T = Float64
p = 1
pl = [[1,2]]
plx = [[1,2]] ######

function dynamics(model::DoubleIntegrator, x, u)
    mp = model.mp  # mass of the point mass in kg (10)
    p = model.p  # number of players
    q1 = x[ @SVector [1,2] ]
    qd1 = x[ @SVector [3,4] ]
    control1 = @SVector [u[pl_ind] for pl_ind in pl[1]]
    qdd1 = control1/mp
    return [qd1; qdd1]
end

# Discretization info
tf = 2.0  # final time
N = 41   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = [SVector{n}([ 0.12, -0.40,  0.00, 0.00]), #lane1
      SVector{n}([-1.00, -0.12,  1.70, 0.00])]
# x0 = [SVector{n}([ 0.12, -0.12,  0.00, 0.00]), #lane1
#       SVector{n}([-1.00, -0.12,  0.00, 0.00])]
xf = [SVector{n}([-0.50,  0.12, -0.30, 0.00]), #lane1
      SVector{n}([ 0.50, -0.12,  1.50, 0.00])]
# xf = [SVector{n}([ 0.12, -0.12,  0.00, 0.00]), #lane1
#       SVector{n}([-0.12, -0.50,  0.00, 0.00])]



dxf = [SVector{n}([1.00, 0.00,  0.00, 0.00]), #lane1
       SVector{n}([1.00, 0.00,  0.00, 0.00])]

# Define a quadratic cost
diag_Q = [SVector{n}([1., 1., 1., 1.]), SVector{n}([10., 10., 1., 0.])]
Q = [0.1*Diagonal(diag_Q[1]),
     0.1*Diagonal(diag_Q[2])]
Qf = [1.0*Diagonal(diag_Q[1]),
      1.0*Diagonal(diag_Q[2])]
R = [0.1*Diagonal(@SVector ones(m)),
  0.1*Diagonal(@SVector ones(m)),
  ]
# Q = [0.1*Diagonal(@SVector ones(n)),
#      0.1*Diagonal(@SVector ones(n)),
#      0.1*Diagonal(@SVector ones(n))]
# Qf = [1.0*Diagonal(@SVector ones(n)),
#      1.0*Diagonal(@SVector ones(n)),
#      1.0*Diagonal(@SVector ones(n))]
# R = [0.1*Diagonal(@SVector ones(length(pl[1]))),
#     0.1*Diagonal(@SVector ones(length(pl[2]))),
#     0.1*Diagonal(@SVector ones(length(pl[3]))),
#     ]
obj = [LQRObjective(Q[i],R[i],Qf[i],xf[i],N) for i=1:2]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = [SVector{m}(zeros(m)) for i=1:2]
U0 = [TO.copy(u0) for k = 1:N]
xs = [SVector{n}(zeros(n)) for i=1:2]
us = [SVector{m}(u0[i]) for i=1:2]
Z = [[KnotPoint(xs[i],us[i],dt) for k = 1:N] for i=1:2]
Z[1][end] = KnotPoint(xs[1],m)
Z[2][end] = KnotPoint(xs[2],m)


# Build problem
car_radius = 0.80/10.
car_radii = [car_radius for i=1:p]
top_road_length = 2.00
bottom_road_length = 0.30
road_width = 0.45
cross_width = 0.15
bound_radius = 0.05
# ramp_length = 1.2
# ramp_angle = pi/12
# obs_radius = 2.0/10. - car_radius
# obs = [0.0, -0.20]
scenario = TIntersectionScenario(
    top_road_length,
    road_width,
    bottom_road_length,
    cross_width,
    car_radii,
    bound_radius)
lanes = [[3], [1]]
# Collision Avoidance
conSet_leader = ConstraintSet(n,m,N)
add_collision_avoidance(conSet_leader, car_radii, plx, p)
add_scenario_constraints(conSet_leader, scenario, lanes[1], plx, N, n; constraint_type=:constraint)
prob_leader = Problem(model, obj[1], conSet_leader, x0[1], xf[1], Z[1], N, tf)
opts_leader = ALTROSolverOptions{T}(projected_newton=false)
altro_leader = ALTROSolver(prob_leader, opts_leader)

@time solve!(altro_leader)

conSet_follower = ConstraintSet(n,m,N)
add_collision_avoidance(conSet_follower, car_radii, plx, p)
add_scenario_constraints(conSet_follower, scenario, lanes[2], plx, N, n; constraint_type=:constraint)
add_leader_constraints(
    conSet_follower,
    TO.states(altro_leader),
    [car_radius, car_radius],
    [[1,2], [1,2]])
prob_follower = Problem(model, obj[2], conSet_follower, x0[2], xf[2], Z[2], N, tf)
opts_follower = ALTROSolverOptions{T}(projected_newton=false)
altro_follower = ALTROSolver(prob_follower, opts_follower)

@time solve!(altro_follower)



include("../src/utils/plot_visualization.jl")
using Plots

X = TO.states(altro_leader)
U = TO.controls(altro_leader)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(altro_leader,pl,plx)

X = TO.states(altro_follower)
U = TO.controls(altro_follower)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(altro_follower,pl,plx)


# animation(altro_leader, scenario, plx)
# animation(altro_follower, scenario, plx)

animation(altro_leader, altro_follower, scenario,
    plx, [car_radius, car_radius]; no_background=true)



a = 100
a = 100
a = 100
a = 100
a = 100


for con in altro_leader.solver_al.solver_uncon.obj.constraints
    @show con.con
    @show con.inds
    @show con.vals
end



for con in altro_follower.solver_al.solver_uncon.obj.constraints
    @show con.con
    @show con.inds
    @show con.vals
end
