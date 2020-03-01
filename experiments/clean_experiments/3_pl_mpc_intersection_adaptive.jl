using ALGAMES
a = 10
# using StaticArrays
# using LinearAlgebra
# const GS = GameSolver


#
# using BenchmarkTools
# using Blink
# using Colors: RGBA, RGB
# using CoordinateTransformations
# using Dates
# using FileIO
# using GeometryTypes
# using JLD2
# using LinearAlgebra
# using Logging
# using MeshCat
# using MeshIO
# using Parameters
# using PartedArrays
# using PGFPlotsX
# using Plots
# using Random
# using SparseArrays
# using StaticArrays
# using Statistics
# using StatsBase
# using Test
# using TrajectoryOptimization
# const TO = TrajectoryOptimization
# const GS = TrajectoryOptimization
#
# using TrajectoryOptimization.Dynamics
# using TrajectoryOptimization.Problems
#
#
# include("../../src/solvers/game_model.jl")
# include("../../src/solvers/game_problem.jl")
# include("../../src/solvers/cost_helpers.jl")
#
# include("../../src/solvers/direct/direct_helpers.jl")
#
# include("../../src/solvers/direct/direct_solver.jl")
# include("../../src/solvers/direct/direct_methods.jl")
# include("../../src/solvers/direct/direct_core.jl")
# include("../../src/solvers/direct/newton_gradient.jl")
# include("../../src/solvers/direct/newton_hessian.jl")
# include("../../src/solvers/inds_helpers.jl")
#
#
# include("../../src/solvers/riccati/algames/algames_solver.jl")
# include("../../src/solvers/riccati/algames/algames_methods.jl")
# include("../../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
# include("../../src/solvers/riccati/ilqgames/ilqgames_methods.jl")
# include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
# include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")
#
# include("../../src/sampler/monte_carlo_sampler.jl")
# include("../../src/sampler/monte_carlo_methods.jl")
#
# include("../../src/scenarios/scenario.jl")
# include("../../src/scenarios/examples/merging.jl")
# include("../../src/scenarios/examples/straight.jl")
# include("../../src/scenarios/examples/t_intersection.jl")
#
# include("../../src/solvers/MPC/mpc_solver.jl")
# include("../../src/solvers/MPC/mpc_methods.jl")
#
# include("../../src/scenarios/scenario_visualization.jl")
#
# include("../../src/utils/constraints.jl")
# include("../../src/utils/monte_carlo_visualization_latex.jl")
# include("../../src/utils/monte_carlo_visualization.jl")
# include("../../src/utils/plot_visualization.jl")
# include("../../src/utils/tests.jl")
# include("../../src/utils/timing.jl")


# Define the dynamics model of the game.
struct InertialUnicycleGame{T} <: AbstractGameModel
    n::Int
    m::Int
    mp::T
	pu::Vector{Vector{Int}}
	px::Vector{Vector{Int}}
    p::Int
end
InertialUnicycleGame() = InertialUnicycleGame(
    12, 6, 1.0, [[1,2],[3,4],[5,6]], [[1,2],[5,6],[9,10]], 3)
Base.size(::InertialUnicycleGame) = 12,6,[[1,2],[3,4],[5,6]],3

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
               0.43, -0.30,  pi/2, 0.12, #lane4
                ]
xf = @SVector [
                1.30, -0.15,  0.00, 0.60, #lane1
               -0.30,  0.15,  pi, 0.60, #lane2
               0.43,  0.35,  pi/2, 0.30, #lane4
               ]
dxf = @SVector [
               0.60,  0.00,  0.00,  0.00, #lane1
              -0.60,  0.00,  0.00,  0.00, #lane2
               0.00,  0.12,  0.00,  0.00, #lane2
              ]

# Define a quadratic cost
diag_Q1 = @SVector [
    0., 1., 1., 1.,
    0., 0., 0., 0.,
    0., 0., 0., 0.]
diag_Q2 = @SVector [
    0., 0., 0., 0.,
    0., 1., 1., 1.,
    0., 0., 0., 0.]
diag_Q3 = @SVector [
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q1),
     0.1*Diagonal(diag_Q2),
     0.1*Diagonal(diag_Q3)]
Qf = [1.0*Diagonal(diag_Q1),
      1.0*Diagonal(diag_Q2),
      1.0*Diagonal(diag_Q3)]

R = [0.1*Diagonal(@SVector ones(length(pu[1]))),
    0.1*Diagonal(@SVector ones(length(pu[2]))),
    0.1*Diagonal(@SVector ones(length(pu[3]))),
    ]
obj = [GS.LQRObjective(Q[i],R[i],Qf[i],xf,N) for i=1:p]

# Define the initial trajectory
# set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector zeros(m)
U0 = [copy(u0) for k = 1:N]
xs = @SVector zeros(n)
us = SVector{m}(u0)
Z = [GS.KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = GS.KnotPoint(xs,m)

dynamics(model, x0, u0)
test_allocation(dynamics, (model, x0, u0))


# Build problem
pr = 3
car_radius_prog = [0.08 + min(k-1, 3)* 0.01 for k = 1:N]
car_radius = 0.08
car_radii = [car_radius for i=1:p]
car_radii_prog = [[car_radius + (j-1)*0.01 for i=1:p] for j=1:4]
con_inds_prog = [j:j for j=1:pr]
con_inds_prog[end] = pr:N
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
conSet_direct = GS.ConstraintSet(n,m,N)
conSet_penalty = GS.ConstraintSet(n,m,N)
con_inds = 2:N
# add_collision_avoidance(conSet_direct, car_radii, px, p, con_inds)
# add_collision_avoidance(conSet_penalty, car_radii, px, p, con_inds)
for j = 1:pr
	add_collision_avoidance(conSet_direct, car_radii_prog[j], px, p, con_inds_prog[j])
	add_collision_avoidance(conSet_penalty, car_radii_prog[j], px, p, con_inds_prog[j])
end
add_scenario_constraints(conSet_direct, scenario, lanes, px, con_inds; constraint_type=:constraint)
add_scenario_constraints(conSet_penalty, scenario, lanes, px, con_inds; constraint_type=:constraint)

# con = BoundConstraint(n, m; u_min=-0.30)
# add_constraint!(conSet_direct, con, 10:N)
# add_constraint!(conSet_penalty, con, 10:N)




# include("../../src/solvers/direct/direct_solver.jl")
prob_direct = GameProblem(model, obj, conSet_direct, x0, xf, Z, N, tf)
prob_penalty = GameProblem(model, obj, conSet_penalty, x0, xf, Z, N, tf)
opts_directgames = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    info_pattern=:open_loop,
	optimality_constraint_tolerance=1e-2)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
opts_penalty = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    iterations_linesearch=10,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05)
solver_penalty = PenaltyiLQGamesSolver(prob_penalty, opts_penalty)
pen = ones(length(solver_penalty.constraints))*100.0
set_penalty!(solver_penalty, pen)

reset!(solver_penalty, reset_type=:full)
@time GS.solve!(solver_penalty)
reset!(solver_directgames, reset_type=:full)
@time GS.solve!(solver_directgames)


a = 10

@btime timing_solve(solver_directgames)
@btime timing_solve(solver_penalty)
# reset!(solver_directgames, reset_type=:full)
# @profiler GS.solve!(solver_directgames)

a = 10


# rollout!(solver_directgames)
# vis, anim = still_animation(solver_directgames, scenario, plx)
# # @time solve!(solver_directgames)
# # full_reset!(solver_directgames)
# # @time solve!(solver_directgames)
# # full_reset!(solver_directgames)
# # # @profiler solve!(solver_directgames)
# # full_reset!(solver_directgames)
# # @btime timing_solve(solver_directgames)
# vis, anim = GS.plan_animation(solver_directgames, scenario, plx;
# 	no_background=false, display_actors=true)
#
# include("../../src/utils/plot_visualization.jl")
# using Plots

X = TO.states(solver_penalty)
U = TO.controls(solver_penalty)
visualize_state(X)
visualize_control(U,pu)
visualize_trajectory_car(solver_penalty)
visualize_collision_avoidance(solver_penalty)
visualize_circle_collision(solver_penalty)
visualize_boundary_collision(solver_penalty)
# visualize_dynamics(solver_penalty)
# visualize_optimality_merit(solver_penalty)
# visualize_H_cond(solver_penalty)
visualize_α(solver_penalty)
visualize_cmax(solver_penalty)
#
X = TO.states(solver_directgames)
U = TO.controls(solver_directgames)
# visualize_state(X_sam)
# visualize_state(X_mpc[1:154])
visualize_state(X)
visualize_control(U,pu)
visualize_trajectory_car(solver_directgames)
visualize_collision_avoidance(solver_directgames)
visualize_circle_collision(solver_directgames)
visualize_boundary_collision(solver_directgames)
visualize_dynamics(solver_directgames)
visualize_optimality_merit(solver_directgames)
visualize_H_cond(solver_directgames)
visualize_α(solver_directgames)
visualize_cmax(solver_directgames)
#

# vis, anim = animation(solver_directgames, scenario;
# 	# vis, anim,
# 	# open_vis=true,
# 	# display_actors=true,
# 	# display_trajectory=false,
# 	no_background=false)
#
# vis, anim = animation(solver_directgames, scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=false,
# 	display_trajectory=false,
# 	no_background=false)

#
# vis.core.tree.children
# vis.core.tree.children["meshcat"]
# vis.core.tree.children["meshcat"].children["roadway"]
# haskey(vis.core.tree.children["meshcat"].children, "roadway")
# delete!(vis.core.tree.children["meshcat"].children, "roadway")
#
#
# haskey(vis.core.tree.children, "meshcat/roadway")
# delete!(vis.core.tree, "meshcat/roadway")
# vis.core.tree.children
#
#
# # MeshCat.update_tree!(vis)
# a = 10
# a = 10

include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

opts_directgames = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    info_pattern=:open_loop,
	constraint_tolerance=1e-3,
	optimality_constraint_tolerance=1e-2)
solver_directgames = DirectGamesSolver(prob_direct, opts_directgames)
@time solve!(solver_directgames)

reset!(solver_directgames, reset_type=:full)

state_noise = @SVector [
    0.008, 0.008, 2*pi/72, 0.03, #+-10cm, +-10cm, +-5deg, +-2.5%
    0.008, 0.008, 2*pi/72, 0.03,
    0.008, 0.008, 2*pi/72, 0.03]
state_noise *= 5.0
selfish_inds = [9,10,11,12]
selfish_dx = [0., 0.12, 0., 0.]

opts_mpc_gamessolver = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=1000,
	N_mpc=50,
	mpc_tf=4.0,
	min_δt=0.001,
	selfish_inds=selfish_inds,
	selfish_dx=selfish_dx,
	noise=state_noise)
mpc_solver = MPCGamesSolver(solver_directgames, dxf, opts_mpc_gamessolver)
reset!(mpc_solver)
solve!(mpc_solver; wait=false)
mpc_solver.stats.iterations
resample!(mpc_solver)

vis, anim = animation(mpc_solver, scenario;
	# vis, anim,
	# open_vis=true,
	display_actors=true,
	display_trajectory=true,
	no_background=false)

vis, anim = animation(mpc_solver, scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true,
	no_background=false)

mean(mpc_solver.stats.solve_time)
Hz = mpc_solver.stats.iterations/mpc_solver.stats.time
mean(mpc_solver.stats.solve_time)
std(mpc_solver.stats.solve_time)
var(mpc_solver.stats.solve_time)
