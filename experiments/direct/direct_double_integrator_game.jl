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
include("lq_game_visualization.jl")

include("direct/stiff_solver.jl")
include("direct/stiff_methods.jl")
include("direct/stiff_helpers.jl")
include("direct/stiff_core.jl")

# Define the dynamics model of the game.
struct DoubleIntegratorGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Array{Array{Int,1},1}
    p::Int
end

DoubleIntegratorGame() = DoubleIntegratorGame(4, 4, 1.0, [[1,2], [3,4]], 2)
Base.size(::DoubleIntegratorGame) = 4,4,[[1,2], [3,4]],2

# Instantiate dynamics model
model = DoubleIntegratorGame()
n,m,pl,p = size(model)

function dynamics(model::DoubleIntegratorGame, x, u)
    mp = model.mp  # mass of the point mass in kg (10)
    p = model.p  # number of players
    pl = model.pl  # control vector partition for each player
    q = x[ @SVector [1,2] ]
    qd = x[ @SVector [3,4] ]
    control = @SVector [u[pl[1][i]] + u[pl[2][i]] for i=1:length(pl[1])]
    qdd = control/mp
    return [qd; qdd]
end


# Discretization info
tf = 2.0  # final time
N = 21   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [1., 2., 3., 4.]
xf = @SVector [0., 0., 0., 0.] # (ie, swing up)

# Define a quadratic cost
Q = [10.0*Diagonal(@SVector ones(n)), -1.0*Diagonal(@SVector ones(n))]
Qf = [100.0*Diagonal(@SVector ones(n)), -10.0*Diagonal(@SVector ones(n))]
R = 10.0*Diagonal(@SVector ones(m))
obj = [LQRObjective(Q[i],R[pl[i],pl[i]],Qf[i],xf,N) for i=1:p]


# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector [0.0, 0.0, 0.0, 0.0]
U0 = [TO.copy(u0) for k = 1:N]
xs = NaN*@SVector zeros(n)
us = SVector{m}(u0)
Zz = [KnotPoint(xs,us,dt) for k = 1:N]
Zz[end] = KnotPoint(xs,m)

# Build problem
T = Float64

conSet_direct = ConstraintSet(n,m,N)
x_min = -100.0* @SVector ones(n)
x_max =  100.0* @SVector ones(n)
u_min = -100.0* @SVector ones(m)
u_max =  100.0* @SVector ones(m)
bound_con = BoundConstraint(n, m, x_min=x_min, x_max=x_max, u_min=u_min, u_max=u_max)
add_constraint!(conSet_direct, bound_con, 1:N)

prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Zz, N, tf)
opts_directgames = DirectGamesSolverOptions{Float64}(
    iterations=10, info_pattern=:open_loop)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)

initial_controls!(prob_direct,U0)
@time solve!(solver_directgames)


cost_expansion(solver_directgames.C, solver_directgames.obj, solver_directgames.Z,
    solver_directgames.model.pl, solver_directgames.model.p)
TO.discrete_jacobian!(solver_directgames.∇F, solver_directgames.model, solver_directgames.Z)
TO.update_active_set!(solver_directgames.constraints,solver_directgames.Z)
TO.evaluate!(solver_directgames.constraints, solver_directgames.Z)
TO.jacobian!(solver_directgames.constraints, solver_directgames.Z)
TO.update_active_set!(solver_directgames.dyn_constraints,solver_directgames.Z) #### maybe not useful
TO.evaluate!(solver_directgames.dyn_constraints, solver_directgames.Z)
update_g_!(solver_directgames)
@show norm(solver_directgames.g_)



# visualize_trajectory(X,U/10,pl)
# visualize_trajectory(X_d,U_d/10,pl)
# visualize_state(X)
# visualize_state(X_d)
# visualize_control(U,pl)
# visualize_control(U_d,pl)

for k=1:N
    println("δ = ", maximum(solver_directgames.Z[k].z - solver_ilqgames.Z[k].z))
end

a = 100
a = 100
a = 100
a = 100
a = 100
