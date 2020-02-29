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

include("ilqgames/lq_game_model.jl")
include("ilqgames/lq_game_problem.jl")
include("ilqgames/lq_game_solver.jl")
include("ilqgames/lq_game_methods.jl")
include("ilqgames/lq_game_helpers.jl")
# include("lq_game_tests.jl")
# include("lq_game_visualization.jl")

include("direct/stiff_solver.jl")
include("direct/stiff_methods.jl")
include("direct/stiff_helpers.jl")
include("direct/stiff_core.jl")
include("direct/stiff_constraints.jl")

# Define the dynamics model of the game.
struct DoubleIntegratorGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Vector{Vector{Int}}
    p::Int
end
# pl1 = @SVector [1,2]
# pl2 = @SVector [3,4]
# pl = [pl1, pl2]
DoubleIntegratorGame() = DoubleIntegratorGame(12, 6, 1.0, [[1,2],[3,4],[5,6]], 3)
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
tf = 4.0  # final time
N = 21   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [-0.5,  -0.05, 0.0, 0.0,
                0.5,   0.05, 0.0, 0.0,
                0.05, -0.5,  0.0, 0.0]
xf = @SVector [ 0.5,  -0.05, 0.0, 0.0,
                -0.05,-0.50, 0.0, 0.0,
                -0.5 , 0.05, 0.0, 0.0]

# Define a quadratic cost
Q = [0.1*Diagonal(@SVector ones(n)),
     0.1*Diagonal(@SVector ones(n)),
     0.1*Diagonal(@SVector ones(n))]
Qf = [1.0*Diagonal(@SVector ones(n)),
      1.0*Diagonal(@SVector ones(n)),
      1.0*Diagonal(@SVector ones(n))]
R = 0.1*Diagonal(@SVector ones(m))
obj = [LQRObjective(Q[i],R[pl[i],pl[i]],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector [0., 0., 0., 0., 0., 0.]
U0 = [TO.copy(u0) for k = 1:N]
xs = NaN*@SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
conSet_direct = ConstraintSet(n,m,N)
# Collision Avoidance
car_radius = 0.05
radius1 = @SVector fill(car_radius, 1)
radius2 = @SVector fill(car_radius, 1)
radius3 = @SVector fill(car_radius, 1)
col_con1 = CollisionConstraint(n, radius1, radius2, 1, 2, 5, 6)
col_con2 = CollisionConstraint(n, radius1, radius3, 1, 2, 9, 10)
col_con3 = CollisionConstraint(n, radius2, radius3, 5, 6, 9, 10)
add_constraint!(conSet_direct, col_con1, 1:N)
add_constraint!(conSet_direct, col_con2, 1:N)
add_constraint!(conSet_direct, col_con3, 1:N)

# Road boundaries
rad = 0.1
B = [SVector{2}([-2., 6.]),
     SVector{2}([ 2., 6.]),
     SVector{2}([-6., 2.]),
     SVector{2}([-2., 2.]),
     SVector{2}([ 2., 2.]),
     SVector{2}([ 6., 2.]),
     SVector{2}([-6.,-2.]),
     SVector{2}([-2.,-2.]),
     SVector{2}([ 2.,-2.]),
     SVector{2}([ 6.,-2.]),
     SVector{2}([-2.,-6.]),
     SVector{2}([ 2.,-6.])]/10.
V = [SVector{2}([1., 0.]),
     SVector{2}([0., 1.])]
M = 10
for inds in plx
    add_circle_boundary!(conSet_direct,inds,Array(B[2]), Array(B[5]), rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[5]), Array(B[6]), rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[10]),Array(B[9]), rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[9]), Array(B[12]),rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[11]),Array(B[8]), rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[8]), Array(B[7]), rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[3]), Array(B[4]), rad,M,n,N)
    add_circle_boundary!(conSet_direct,inds,Array(B[4]), Array(B[1]), rad,M,n,N)
end

for inds in plx
    con1 = BoundaryConstraint(n, B[2],  B[5],  V[1], inds[1], inds[2])
    con2 = BoundaryConstraint(n, B[5],  B[6],  V[2], inds[1], inds[2])
    con3 = BoundaryConstraint(n, B[10], B[9], -V[2], inds[1], inds[2])
    con4 = BoundaryConstraint(n, B[9],  B[12], V[1], inds[1], inds[2])
    con5 = BoundaryConstraint(n, B[11], B[8], -V[1], inds[1], inds[2])
    con6 = BoundaryConstraint(n, B[8],  B[7], -V[2], inds[1], inds[2])
    con7 = BoundaryConstraint(n, B[3],  B[4],  V[2], inds[1], inds[2])
    con8 = BoundaryConstraint(n, B[4],  B[1], -V[1], inds[1], inds[2])
    add_constraint!(conSet_direct, con1, 1:N)
    add_constraint!(conSet_direct, con2, 1:N)
    add_constraint!(conSet_direct, con3, 1:N)
    add_constraint!(conSet_direct, con4, 1:N)
    add_constraint!(conSet_direct, con5, 1:N)
    add_constraint!(conSet_direct, con6, 1:N)
    add_constraint!(conSet_direct, con7, 1:N)
    add_constraint!(conSet_direct, con8, 1:N)
end


include("direct/stiff_solver.jl")
include("direct/stiff_methods.jl")
include("direct/stiff_helpers.jl")
include("direct/stiff_core.jl")
include("direct/stiff_constraints.jl")

prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{Float64}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=5,
    info_pattern=:open_loop)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)

initial_controls!(prob_direct,U0)

# rollout!(solver_directgames)
# cost_expansion(solver_directgames.C, solver_directgames.obj, solver_directgames.Z,
#     solver_directgames.model.pl, solver_directgames.model.p)
# TO.update_active_set!(solver_directgames.constraints,solver_directgames.Z)
# TO.evaluate!(solver_directgames.constraints, solver_directgames.Z)
# TO.jacobian!(solver_directgames.constraints, solver_directgames.Z)
# TO.discrete_jacobian!(solver_directgames.∇F, solver_directgames.model,
#     solver_directgames.Z)
# TO.evaluate!(solver_directgames.dyn_constraints, solver_directgames.Z)
#
# update_H_!(solver_directgames)

# solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
# @time solve!(solver_directgames)
@time solve!(solver_directgames)
reset!(solver_directgames, prob_direct)
@time solve!(solver_directgames)


function timing_solve(solver::DirectGamesSolver{T}, prob::GameProblem) where T
    reset!(solver_directgames, prob_direct)
    solve!(solver_directgames)
end
function timing_reset(solver::DirectGamesSolver{T}, prob::GameProblem) where T
    reset!(solver_directgames, prob_direct)
end

@btime timing_solve(solver_directgames, prob_direct)
@btime timing_reset(solver_directgames, prob_direct)


a = 100
a = 100


# cost_expansion(solver_directgames.C, solver_directgames.obj, solver_directgames.Z,
#     solver_directgames.model.pl, solver_directgames.model.p)
# TO.update_active_set!(solver_directgames.constraints,solver_directgames.Z)
# TO.evaluate!(solver_directgames.constraints, solver_directgames.Z)
# TO.jacobian!(solver_directgames.constraints, solver_directgames.Z)
# TO.discrete_jacobian!(solver_directgames.∇F, solver_directgames.model,
#     solver_directgames.Z)
# TO.evaluate!(solver_directgames.dyn_constraints, solver_directgames.Z)
# update_g_!(solver_directgames)
# @show norm(solver_directgames.g_)

include("lq_game_visualization.jl")
using Plots
X = TO.states(solver_directgames)
U = TO.controls(solver_directgames)
visualize_state(X)
visualize_control(U,pl)
visualize_trajectory_car(solver_directgames,X,U/5,pl,plx,car_radius)
visualize_collision(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)
# for k=1:N
#     println("δ = ", maximum(solver_directgames.Z[k].z - solver_ilqgames.Z[k].z))
# end

a = 100
a = 100
a = 100
a = 100
a = 100




#
# T = Float64
# st = DirectGamesStats{T}()
# reset!(st, 10, 20, 2)
# record_iteration!(st, [100.,200.], [10.,20.], 0.1)
# st
# st.cost
# st.dJ
# st.cmax
# record_iteration!(st, [101.,201.], [11.,21.], 0.11)
#

#
# pls = []
# for i in 1:length(a_pl)
#     @show i
#     a = SVector{length(a_pl[i])}(a_pl[i])
#     push!(pls, a)
# end
# pls
# a = SVector{length(a_pl[1])}(a_pl[1])
#
#
# temp = []
# for inds in a_pl
#     push!(temp, SVector{length(inds)}(inds))
# end
# pls = [inds for inds in temp]
#
