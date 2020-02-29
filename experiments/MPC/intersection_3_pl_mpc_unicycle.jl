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

include("../src/solvers/MPC/mpc_solver.jl")
include("../src/solvers/MPC/mpc_methods.jl")

include("../src/utils/tests.jl")

# Define the dynamics model of the game.
struct UnicycleGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
    pl::Vector{Vector{Int}}
    p::Int
end
UnicycleGame() = UnicycleGame(
    9, 6, 1.0, [[1,2],[3,4],[5,6]], 3)
Base.size(::UnicycleGame) = 9,6,[[1,2],[3,4],[5,6]],3

# Instantiate dynamics model
model = UnicycleGame()
n,m,pl,p = size(model)
T = Float64
plx = [[1,2], [4,5], [7,8]]


function dynamics(model::UnicycleGame, x, u)
    qd1 = @SVector [cos(x[3]*2*pi), sin(x[3]*2*pi)]
    qd1 *= u[1]
    qd2 = @SVector [cos(x[6]*2*pi), sin(x[6]*2*pi)]
    qd2 *= u[3]
    qd3 = @SVector [cos(x[9]*2*pi), sin(x[9]*2*pi)]
    qd3 *= u[5]
    ω1 = u[ @SVector [2] ]
    ω2 = u[ @SVector [4] ]
    ω3 = u[ @SVector [6] ]
    return [qd1; ω1; qd2; ω2; qd3; ω3]
end

# Discretization info
tf = 2.0  # final time
N = 21   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
# x0 = @SVector [
#                -0.15,  0.3,  -pi/2, #lane1
#                 0.4,   0.15,  pi, #lane2
#                -0.9,  -0.15,  0.0, #lane4
#                 ]
# xf = @SVector [
#                 0.6,  -0.15, 0.0,    #lane1
#                -0.15, -0.50, 3*pi/2, #lane2
#                # -0.15 ,  0.3,  0.0, #lane4
#                 0.15,  0.20, pi/2,   #lane4
#                ]

x0 = @SVector [
              -0.15,  0.3,  -0.25, #lane1
               0.6,   0.15,  0.5, #lane2
              -1.0,  -0.15,  0.0, #lane4
               ]
xf = @SVector [
               0.6,  -0.15,  0.0,    #lane1
               # -0.15, -0.50, 0.75, #lane2
              -0.15, -0.20,  0.75, #lane2
              # 0.15,  0.20, 0.25,   #lane4
               0.15,  0.50,  0.25,   #lane4
              ]
uf = @SVector [
               [0.5, 0.0], #lane1
               [0.5, 0.0], #lane2
               [0.5, 0.0],   #lane4
               ]

# Define a quadratic cost
Q = [0.1*Diagonal(@SVector ones(n)),
     0.1*Diagonal(@SVector ones(n)),
     0.1*Diagonal(@SVector ones(n))]
Qf = [1.0*Diagonal(@SVector ones(n)),
      1.0*Diagonal(@SVector ones(n)),
      1.0*Diagonal(@SVector ones(n))]
R = 0.1*Diagonal(@SVector ones(m))
obj = [LQRObjective(Q[i],R[pl[i],pl[i]],Qf[i],xf,uf[i],N) for i=1:p]





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
car_radius = 0.8/10. # 0.8/10.
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
cir_rad = 0.5/10.
l1 = 3.0/10.
l2 = 7.0/10.

add_intersection_constraints(conSet_direct, l1, l2, cir_rad, car_radius,
    1, plx, N, lane=1)
add_intersection_constraints(conSet_direct, l1, l2, cir_rad, car_radius,
    2, plx, N, lane=2)
add_intersection_constraints(conSet_direct, l1, l2, cir_rad, car_radius,
    3, plx, N, lane=4)

include("../src/solvers/direct/direct_solver.jl")
include("../src/solvers/direct/direct_methods.jl")
include("../src/solvers/direct/direct_helpers.jl")
include("../src/solvers/direct/direct_core.jl")
include("../src/solvers/direct/direct_constraints.jl")

prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{Float64}(
    iterations=10,
    inner_iterations=15,
    iterations_linesearch=10,
    # compute_cond=true,
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


a = 100
a = 100
a = 100
a = 100
a = 100







include("../src/solvers/direct/direct_solver.jl")
include("../src/solvers/direct/direct_methods.jl")


include("../src/solvers/MPC/mpc_solver.jl")
include("../src/solvers/MPC/mpc_methods.jl")


solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
initial_controls!(prob_direct,U0)
@time solve!(solver_directgames)
full_reset!(solver_directgames, prob_direct)
mpc_solver = MPCGamesSolver(solver_directgames)
reset!(mpc_solver)
solve!(mpc_solver)



discrete_dynamics(
    TO.DEFAULT_Q,
    mpc_solver.solver.model,
    TO.state(mpc_solver.solver.Z[1]),
    TO.control(mpc_solver.solver.Z[1]),
    0.0,
    0.1)



# i = TO.Explicit
# L = n+m
#
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


for con in solver_directgames.dyn_constraints.constraints
    # @show typeof(con.con)
    @show typeof(con.con)
    @show maximum(con.c_max)
end
for con in solver_directgames.constraints.constraints
    # @show typeof(con.con)
    # @show typeof(con.con)
    @show maximum(con.c_max)
end


@show (solver_directgames.dyn_constraints.constraints[1].vals[end])
for

a = 100
a = 100


#
# function generate_intersection(i1, i2, circle_radius, car_radius)
#     B = [SVector{2}([-i1, i2]),
#          SVector{2}([ i1, i2]),
#          SVector{2}([-i2, i1]),
#          SVector{2}([-i1, i1]),
#          SVector{2}([ i1, i1]),
#          SVector{2}([ i2, i1]),
#          SVector{2}([-i2,-i1]),
#          SVector{2}([-i1,-i1]),
#          SVector{2}([ i1,-i1]),
#          SVector{2}([ i2,-i1]),
#          SVector{2}([-i1,-i2]),
#          SVector{2}([ i1,-i2])]
#     i_ = i1 + circle_radius - car_radius
#     C = [SVector{2}([-i_, i2]),
#          SVector{2}([ i_, i2]),
#          SVector{2}([-i2, i_]),
#          SVector{2}([-i_, i_]),
#          SVector{2}([ i_, i_]),
#          SVector{2}([ i2, i_]),
#          SVector{2}([-i2,-i_]),
#          SVector{2}([-i_,-i_]),
#          SVector{2}([ i_,-i_]),
#          SVector{2}([ i2,-i_]),
#          SVector{2}([-i_,-i2]),
#          SVector{2}([ i_,-i2])]
#     return B,C
# end
# B, C = generate_intersection(i1, i2, cir_rad, car_radius)
# # B = [SVector{2}([-1.5, 6.0]),
# #      SVector{2}([ 1.5, 6.0]),
# #      SVector{2}([-6.0, 1.5]),
# #      SVector{2}([-1.5, 1.5]),
# #      SVector{2}([ 1.5, 1.5]),
# #      SVector{2}([ 6.0, 1.5]),
# #      SVector{2}([-6.0,-1.5]),
# #      SVector{2}([-1.5,-1.5]),
# #      SVector{2}([ 1.5,-1.5]),
# #      SVector{2}([ 6.0,-1.5]),
# #      SVector{2}([-1.5,-6.0]),
# #      SVector{2}([ 1.5,-6.0])]/10.
# # C = [SVector{2}([-1.5 - (rad-car_radius) , 6.0]),
# #      SVector{2}([ 1.5, 6.0]),
# #      SVector{2}([-6.0, 1.5]),
# #      SVector{2}([-1.5, 1.5]),
# #      SVector{2}([ 1.5, 1.5]),
# #      SVector{2}([ 6.0, 1.5]),
# #      SVector{2}([-6.0,-1.5]),
# #      SVector{2}([-1.5,-1.5]),
# #      SVector{2}([ 1.5,-1.5]),
# #      SVector{2}([ 6.0,-1.5]),
# #      SVector{2}([-1.5,-6.0]),
# #      SVector{2}([ 1.5,-6.0])]/10.
#
# V = [SVector{2}([1., 0.]),
#      SVector{2}([0., 1.])]
# M = 10
# for inds in plx
#     add_circle_boundary!(conSet_direct,inds,Array(C[2]), Array(C[5]), cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[5]), Array(C[6]), cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[10]),Array(C[9]), cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[9]), Array(C[12]),cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[11]),Array(C[8]), cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[8]), Array(C[7]), cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[3]), Array(C[4]), cir_rad,M,n,N)
#     add_circle_boundary!(conSet_direct,inds,Array(C[4]), Array(C[1]), cir_rad,M,n,N)
# end
#
# for inds in plx
#     con1 = BoundaryConstraint(n, B[2],  B[5],  V[1], inds[1], inds[2])
#     con2 = BoundaryConstraint(n, B[5],  B[6],  V[2], inds[1], inds[2])
#     con3 = BoundaryConstraint(n, B[10], B[9], -V[2], inds[1], inds[2])
#     con4 = BoundaryConstraint(n, B[9],  B[12], V[1], inds[1], inds[2])
#     con5 = BoundaryConstraint(n, B[11], B[8], -V[1], inds[1], inds[2])
#     con6 = BoundaryConstraint(n, B[8],  B[7], -V[2], inds[1], inds[2])
#     con7 = BoundaryConstraint(n, B[3],  B[4],  V[2], inds[1], inds[2])
#     con8 = BoundaryConstraint(n, B[4],  B[1], -V[1], inds[1], inds[2])
#     add_constraint!(conSet_direct, con1, 1:N)
#     add_constraint!(conSet_direct, con2, 1:N)
#     add_constraint!(conSet_direct, con3, 1:N)
#     add_constraint!(conSet_direct, con4, 1:N)
#     add_constraint!(conSet_direct, con5, 1:N)
#     add_constraint!(conSet_direct, con6, 1:N)
#     add_constraint!(conSet_direct, con7, 1:N)
#     add_constraint!(conSet_direct, con8, 1:N)
# end


k = 2
i = 1
con = solver_directgames.constraints.constraints[end]
rel_zind_i = rel_zinds(con,solver_directgames.Z[k],k,i,n,m,pl,p,N)
zind_i = zinds(solver_directgames,con,k,i)
∇ci = con.∇c[k][:,rel_zind_i]
Iμ = TO.penalty_matrix(con,k)
@btime $solver_directgames.g_[$zind_i] += $∇ci'*$con.λ[$k] + $∇ci'*$Iμ*$con.vals[$k]
@btime $solver_directgames.state_view_g[1][2] += $∇ci'*$con.λ[$k] + $∇ci'*$Iμ*$con.vals[$k]
@btime $∇ci'*$con.λ[$k] + $∇ci'*$Iμ*$con.vals[$k]

@show 1


@btime rel_zind_i = rel_zinds(con,solver_directgames.Z[k],k,i,n,m,pl,p,N)
rel_zind = rel_zinds(con,solver_directgames.Z[k],k,n,m,pl,p,N)
@btime zind_i = zinds(solver_directgames,con,k,i)
zind = zinds(solver_directgames,con,k)
@btime ∇ci = con.∇c[k][:,rel_zind_i] #ok
∇c = con.∇c[k][:,rel_zind] #ok
@btime Iμ = TO.penalty_matrix(con,k)
@btime solver_directgames.H_[zind_i,zind] .= ∇ci'*Iμ*∇c
@btime ∇ci'*Iμ*∇c
@btime solver_directgames.H_[zind_i,zind]
vH_ = view(H, zind_i, zind)
@btime vH_ .+= ∇ci'*Iμ*∇c
@btime solver_directgames.state_view_H[1][2] .+= ∇ci'*Iμ*∇c
H = Array(solver_directgames.H_)
@btime H[zind_i, zind]


@btime $rel_zind_i = rel_zinds($con,$solver_directgames.Z[$k],$k,$i,$n,$m,$pl,$p,$N)
@btime rel_zind = rel_zinds($con,$solver_directgames.Z[k],$k,$n,$m,$pl,$p,$N)
@btime $zind_i = zinds($solver_directgames,$con,$k,$i)
@btime $zind = zinds($solver_directgames,$con,$k)
@btime $∇ci = $con.∇c[$k][:,$rel_zind_i] #ok
@btime $∇c = $con.∇c[$k][:,$rel_zind] #ok
@btime $Iμ = TO.penalty_matrix($con,$k)
@btime $solver_directgames.H_[$zind_i,$zind] .+= $∇ci'*$Iμ*$∇c
@btime $∇ci'*$Iμ*$∇c
@btime solver_directgames.H_[$zind_i,$zind]
vH_ = view(H, zind_i, zind)
@btime $vH_ .+= $∇ci'*$Iμ*$∇c
@btime $solver_directgames.state_view_H[1][2] .+= $∇ci'*$Iμ*$∇c
H = Array(solver_directgames.H_)
@btime $H[$zind_i, $zind]


@btime $solver_directgames.H_[$zind_i,$zind] .+= $∇ci'*$Iμ*$∇c
@btime $solver_directgames.state_view_H[1][2] .+= $∇ci'*$Iμ*$∇c


struct mystruct7{T}
    a::Int
    b::T
    vH::SubArray{T,2,SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}
end

ms = mystruct(1, 2.3)
ms2 = mystruct3(1, 2.3, vH_)
ms2 = mystruct4(1, 2.3, vH_)
ms2 = mystruct5(1, 2.3, vH_)
ms2 = mystruct6(1, 2.3, vH_)
ms2 = mystruct7(1, 2.3, vH_)
