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
include("lq_game_tests.jl")
include("lq_game_visualization.jl")

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

# Instantiate dynamics model
model = DoubleIntegratorGame()
n,m,pl,p = size(model)


# Discretization info
tf = 2.0  # final time
N = 21   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [1., 2., 3., 4.]
xf = @SVector [0., 0., 0., 0.] # (ie, swing up)

# Define a quadratic cost
# Q = 10.0*Diagonal(@SVector ones(n))
# Qf = 100.0*Diagonal(@SVector ones(n))
# R = 1.0*Diagonal(@SVector ones(m))
Q = [10.0*Diagonal(@SVector ones(n)), -1.0*Diagonal(@SVector ones(n))]
Qf = [100.0*Diagonal(@SVector ones(n)), -10.0*Diagonal(@SVector ones(n))]
R = 10.0*Diagonal(@SVector ones(m))
obj = [LQRObjective(Q[i],R[pl[i],pl[i]],Qf[i],xf,N) for i=1:p]

# Set Constraints
conSet = ConstraintSet(n,m,N)

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector [0.0, 0.0, 0.0, 0.0]
U0 = [TO.copy(u0) for k = 1:N]
xs = NaN*@SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
#   TODO: better constructors
T = Float64
prob = GameProblem(model, obj, conSet, x0, xf, Z, N, tf)
# prob = GameProblem{RK3}(model, obj, conSet, x0, xf, Z, N, tf)

# Build the solver
opts_ilqgames = iLQGamesSolverOptions{Float64}(
    iterations=10, info_pattern=:open_loop)
# opts_ilqgames = iLQGamesSolverOptions{Float64}(
    # iterations=10, info_pattern=:feedback)
solver_ilqgames = iLQGamesSolver(prob, opts_ilqgames)

# Convert to Augmented Lagrangian problem
#   TODO: convert this automatically (will incur allocations)
# prob_al = TO.convertProblem(prob, solver)

# Solve
#   reset the control before the solve so you can re-run it easily
initial_controls!(prob, u0)
solve!(solver_ilqgames)


# Analyze results
#   look in solver.stats for stats on the solve
println("Outer loop iterations: ", solver_ilqgames.stats.iterations)

# #   get final cost
# println("Final cost: ", cost(prob))
#   extract trajectories and plot
X = TO.states(solver_ilqgames)
U = TO.controls(solver_ilqgames)
# plot(X, 1:4)
# plot(U)
# visualize_trajectory(X,U/10,pl)
# visualize_state(X)
# visualize_control(U,pl)

# colinear = [U[k][pl[1]]'*U[k][pl[2]] / norm(U[k][pl[1]]) / norm(U[k][pl[2]]) for k=1:N-1]
# plot([[0];colinear])




#
# #   plot cost convergence
# plot(1:solver_ilqgames.stats.iterations, solver_ilqgames.stats.cost, yscale=:log10)
#
# # Benchmark the result
# @btime begin
#     initial_controls!(prob, U0)
#     solve!(prob, solver_ilqgames)
# end
#
# there shouldn't be any allocations
TO.initial_controls!(solver_ilqgames, U0)
@allocated solve!(solver_ilqgames)
@time solve!(solver_ilqgames)
# initial_controls!(prob, U0)
# @btime solve!(prob, solver_ilqgames)







include("direct/stiff_solver.jl")
include("direct/stiff_methods.jl")
include("direct/stiff_helpers.jl")
include("direct/stiff_core.jl")

N_direct = 21
x0_direct = @SVector [1., 2., 3., 4.]
Z_direct = [KnotPoint(xs,us,dt) for k = 1:N_direct]
Z_direct[end] = KnotPoint(xs,m)
obj_direct = [LQRObjective(Q[i],R[pl[i],pl[i]],Qf[i],xf,N_direct) for i=1:p]
conSet_direct = ConstraintSet(n,m,N)
x_min = -100.0* @SVector ones(n)
x_max =  100.0* @SVector ones(n)
u_min = -100.0* @SVector ones(m)
u_max =  100.0* @SVector ones(m)
bound_con = BoundConstraint(n, m, x_min=x_min, x_max=x_max, u_min=u_min, u_max=u_max)
add_constraint!(conSet_direct, bound_con, 1:N)


prob_direct = GameProblem(model, obj_direct, conSet_direct, x0_direct, xf, Z_direct, N_direct, tf)
opts_directgames = DirectGamesSolverOptions{Float64}(
    iterations=10, info_pattern=:open_loop)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)

initial_controls!(prob_direct,U0)
@time solve!(solver_directgames)


cost_expansion(solver_directgames.C, solver_directgames.obj, solver_directgames.Z,
    solver_directgames.model.pl, solver_directgames.model.p)
TO.discrete_jacobian!(solver_directgames.∇F, solver_directgames.model, solver_directgames.Z)
TO.update_active_set!(solver_directgames.constraints,Z)
TO.evaluate!(solver_directgames.constraints, Z)
TO.jacobian!(solver_directgames.constraints, Z)
TO.update_active_set!(solver_directgames.dyn_constraints,Z) #### maybe not useful
TO.evaluate!(solver_directgames.dyn_constraints, Z)
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

# # Find dyn
# dyn_conSet = solver_directgames.dyn_constraints
# dyn_conSet.constraints[1]
#
#
# #Check
# H_ = solver_directgames.H_
# g_ = solver_directgames.g_
#
# ix = solver_directgames.Z[1]._x
# iu = solver_directgames.Z[1]._u
# A = round.(solver_directgames.∇F[1][ix,ix],digits=4)
# B = round.(solver_directgames.∇F[1][ix,iu],digits=4)
#
# row_off = (n*p+m)*(N_direct-1)
# col_off = (n+m)*(N_direct-1)
#
# #DYNAMICS
# @show A
# @show B
#
# #TOP LEFT
# @show H1 = round.(Array(H_[1:row_off,1:col_off]),digits=4)
# #TOP RIGHT
# @show H2 = round.(Array(H_[1:row_off,col_off+1:end]),digits=4)
# #BOTTOM LEFT
# @show H3 = round.(Array(H_[row_off+1:end,1:col_off]),digits=4)
# #BOTTOM RIGHT
# @show H4 = H_[row_off+1:end, col_off+1:end]
#
#
# #TOP
# @show g1 = round.(g_[1:24],digits=4)
# #BOTTOM
# @show g2 = round.(g_[1+24:end],digits=4)


#
#
#
# conSet_direct0 = ConstraintSet(n,m,N)
# x_min = -100.0* @SVector ones(n)
# x_max =  100.0* @SVector ones(n)
# u_min = -100.0* @SVector ones(m)
# u_max =  100.0* @SVector ones(m)
# bound_con0 = BoundConstraint(n, m, x_min=x_min, x_max=x_max, u_min=u_min, u_max=u_max)
# bound_con1 = BoundConstraint(n, m, x_min=x_min/2, x_max=x_max/2, u_min=u_min/2, u_max=u_max/2)
# add_constraint!(conSet_direct0, bound_con0, 1:N)
# add_constraint!(conSet_direct0, bound_con1, 2:N-1)
#
# conLen = length.(conSet_direct0.constraints)
# idx = 0
# for k = 1:N
#     for (i,con) in enumerate(conSet_direct0.constraints)
#         @show con.inds
#         if k in con.inds
#             j = TO._index(con,k)
#             @show i
#             @show j
#             # @show idx
#             @show conLen[i]
#             # cons[i][j] = idx .+ (1:conLen[i])
#             # idx += conLen[i]
#         end
#     end
# end


# conLen = length.(conSet.constraints)
#
# cons = [[@SVector ones(Int,length(con)) for i in eachindex(con.inds)] for con in conSet.constraints]
#
#
# a = 1
# a = 1
# a = 1
# a = 1
# a = 1

# blk_len = map(con->length(con.∇c[1]), conSet_direct0.constraints)
# blk_len = map(con->length(con.∇c[1]), dyn_conSet.constraints)
#
#
# blk_len = map(con->length(con.∇c[1]), conSet_direct0.constraints)
# con_len = map(con->length(con.∇c), conSet_direct0.constraints)
# con_inds_const = TO.gen_con_inds(conSet_direct0, :by_constraint)
# con_inds_knot = TO.gen_con_inds(conSet_direct0, :by_knotpoint)
#
# # These get filled in with call to constraint_jacobian_structure
# linds = [[@SVector zeros(Int,blk) for i = 1:len] for (blk,len) in zip(blk_len, con_len)]
#
#
# linds = TO.jacobian_linear_inds(prob_direct)


# # Linear indices
# if structure == :by_constraint
#     for (i,con) in enumerate(conSet.constraints)
#         for (j,k) in enumerate(con.inds)
#             inds = idx .+ (1:blk_len[i])
#             linds[i][j] = inds
#             con.∇c[j] = inds
#             idx += blk_len[i]
#         end
#     end
# elseif structure == :by_knotpoint
#     for k = 1:N
#         for (i,con) in enumerate(conSet.constraints)
#             if k in con.inds
#                 inds = idx .+ (1:blk_len[i])
#                 j = _index(con,k)
#                 linds[i][j] = inds
#                 con.∇c[j] = inds
#                 idx += blk_len[i]
#             end
#         end
#     end
# end


# Find ∇dyn

# for conVals in solver_directgames.constraints.constraints
#     for (j,k) in enumerate(conVals)
#         @show j, k
#     end
# end
# for conVals in solver_directgames.dyn_constraints.constraints
#     for (j,k) in enumerate(conVals.inds)
#         @show j, k
#     end
# end

#
# if structure == :by_constraint
#     for (i,con) in enumerate(conSet.constraints)
#         for (j,k) in enumerate(con.inds)
#             cons[i][_index(con,k)] = idx .+ (1:conLen[i])
#             idx += conLen[i]
#         end
#     end
# elseif structure == :by_knotpoint
#     for k = 1:N
#         for (i,con) in enumerate(conSet.constraints)
#             if k in con.inds
#                 j = _index(con,k)
#                 cons[i][j] = idx .+ (1:conLen[i])
#                 idx += conLen[i]
#             end
#         end
#     end
# end
