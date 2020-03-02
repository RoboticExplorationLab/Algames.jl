
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

include("../../src/solvers/game_model.jl")
include("../../src/solvers/game_problem.jl")
include("../../src/solvers/cost_helpers.jl")

include("../../src/solvers/direct/direct_helpers.jl")

include("../../src/solvers/direct/direct_solver.jl")
include("../../src/solvers/direct/direct_methods.jl")
include("../../src/solvers/direct/direct_core.jl")
include("../../src/solvers/direct/newton_gradient.jl")
include("../../src/solvers/direct/newton_hessian.jl")
include("../../src/solvers/inds_helpers.jl")


include("../../src/solvers/riccati/algames/algames_solver.jl")
include("../../src/solvers/riccati/algames/algames_methods.jl")
include("../../src/solvers/riccati/ilqgames/ilqgames_solver.jl")
include("../../src/solvers/riccati/ilqgames/ilqgames_methods.jl")
include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_solver.jl")
include("../../src/solvers/riccati/penalty_ilqgames/penalty_ilqgames_methods.jl")

include("../../src/sampler/monte_carlo_sampler.jl")
include("../../src/sampler/monte_carlo_methods.jl")

include("../../src/scenarios/scenario.jl")
include("../../src/scenarios/examples/merging.jl")
include("../../src/scenarios/examples/straight.jl")
include("../../src/scenarios/examples/t_intersection.jl")

include("../../src/solvers/MPC/mpc_solver.jl")
include("../../src/solvers/MPC/mpc_methods.jl")

include("../../src/scenarios/scenario_visualization.jl")
include("../../src/scenarios/adaptive_plot.jl")

include("../../src/utils/constraints.jl")
include("../../src/utils/monte_carlo_visualization_latex.jl")
include("../../src/utils/monte_carlo_visualization.jl")
include("../../src/utils/plot_visualization.jl")
include("../../src/utils/tests.jl")
include("../../src/utils/timing.jl")


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
	12,
	6,
	1.0,
	[[1,2],[3,4],[5,6]],
	[[1,2],[5,6],[9,10]],
	3)
Base.size(::InertialUnicycleGame) = 12,6,[[1,2],[3,4],[5,6]],3 # n,m,pu,p

# Instantiate dynamics model
model = InertialUnicycleGame()
n,m,pu,p = size(model)
T = Float64
px = model.px
function TO.dynamics(model::InertialUnicycleGame, x, u) # Non memory allocating dynamics
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
N = 21    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.80, -0.05,  0.00, 0.90, #lane1
               -1.00, -0.05,  0.00, 0.90, #lane2
               -0.90, -0.30, pi/12, 0.95, #lane4
                ]
xf = @SVector [
                1.10, -0.05,  0.00, 0.90, #lane1
                0.70, -0.05,  0.00, 0.90, #lane2
                0.90, -0.05,  0.00, 0.90, #lane4
               ]

# Define a quadratic cost
diag_Q1 = @SVector [ # Player 1 state cost
    0., 10., 1., 1.,
    0., 0., 0., 0.,
    0., 0., 0., 0.]
diag_Q2 = @SVector [ # Player 2 state cost
    0., 0., 0., 0.,
    0., 10., 1., 1.,
    0., 0., 0., 0.]
diag_Q3 = @SVector [ # Player 3 state cost
    0., 0., 0., 0.,
    0., 0., 0., 0.,
    0., 10., 1., 1.]
Q = [0.1*Diagonal(diag_Q1), # Players state costs
     0.1*Diagonal(diag_Q2),
     0.1*Diagonal(diag_Q3)]
Qf = [1.0*Diagonal(diag_Q1),
      1.0*Diagonal(diag_Q2),
      1.0*Diagonal(diag_Q3)]

# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[1]))),
     0.1*Diagonal(@SVector ones(length(pu[2]))),
     0.1*Diagonal(@SVector ones(length(pu[3]))),
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
actors_types = [:car for i=1:p]
road_length = 2.20
road_width = 0.30
ramp_length = 1.2
ramp_angle = pi/12
scenario = MergingScenario(road_length, road_width, ramp_length,
	ramp_angle, actors_radii, actors_types)

# Create constraints
conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(conSet, actors_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(conSet, scenario, px, con_inds; constraint_type=:constraint)


prob = GameProblem(model, obj, conSet, x0, xf, Z, N, tf)
opts = DirectGamesSolverOptions{T}(
    iterations=10,
	record_condition=false, ##################
    inner_iterations=20,
    iterations_linesearch=10,
	optimality_constraint_tolerance=1e-2,
	log_level=Logging.Debug)
solver = DirectGamesSolver(prob, opts)
reset!(solver, reset_type=:full)
@time solve!(solver)

# @btime timing_solve(solver)



conSet_penalty = ConstraintSet(n,m,N)

add_collision_avoidance(conSet_penalty, actors_radii, px, p, con_inds)

add_scenario_constraints(conSet_penalty, scenario, px, con_inds; constraint_type=:constraint)

prob_penalty = GameProblem(model, obj, conSet_penalty, x0, xf, Z, N, tf)

opts_penalty = PenaltyiLQGamesSolverOptions{T}(
    iterations=200,
    gradient_norm_tolerance=1e-2,
    cost_tolerance=1e-4,
    iterations_linesearch=5,
    line_search_lower_bound=0.0,
    line_search_upper_bound=0.05)
solver_penalty = PenaltyiLQGamesSolver(prob_penalty, opts_penalty)
pen = ones(length(solver_penalty.constraints))*100.0
set_penalty!(solver_penalty, pen)
@time solve!(solver_penalty)
# reset!(solver_directgames, reset_type=:full)
# @profiler GS.solve!(solver_directgames)

X = TO.states(solver)
U = TO.controls(solver)
visualize_state(X)
visualize_control(U,pu)
visualize_trajectory_car(solver)
visualize_collision_avoidance(solver)
visualize_circle_collision(solver)
visualize_boundary_collision(solver)
visualize_dynamics(solver)
visualize_optimality_merit(solver)
visualize_H_cond(solver)
visualize_α(solver)
visualize_cmax(solver)

X = TO.states(solver_penalty)
U = TO.controls(solver_penalty)
visualize_state(X)
visualize_control(U,pu)
visualize_trajectory_car(solver_penalty)
visualize_collision_avoidance(solver_penalty)
visualize_circle_collision(solver_penalty)
visualize_boundary_collision(solver_penalty)
visualize_α(solver_penalty)
visualize_cmax(solver_penalty)




# vis=Visualizer()
# anim=MeshCat.Animation()
# open(vis)
# vis, anim = animation(solver_directgames, scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=false,
# 	no_background=false)
#
#
# vis, anim = animation(solver_directgames, scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=false,
# 	no_background=false)


state_noise = @SVector [
    0.06, 0.06, 2*pi/72, 0.05,
    0.06, 0.06, 2*pi/72, 0.05,
    0.06, 0.06, 2*pi/72, 0.05]
opts_monte_carlo = MonteCarloSamplerOptions{n,T}(
    noise=state_noise,
    iterations=200)
sampler = MonteCarloSampler(solver, opts_monte_carlo)
sampler_penalty = MonteCarloSampler(solver_penalty, opts_monte_carlo)

monte_carlo_sampling(sampler)
monte_carlo_sampling(sampler_penalty)


using PGFPlotsX
using StatsBase
include("../../src/utils/monte_carlo_visualization.jl")
visualize_H_cond(sampler; save=true)
visualize_solve_time(sampler; save=true)
visualize_cmax(sampler; save=true)
visualize_iterations_total(sampler; save=true)
visualize_optimality_merit(sampler; save=true)

visualize_H_cond(sampler_penalty; save=true)
visualize_solve_time(sampler_penalty; save=true)
visualize_cmax(sampler_penalty; save=true)
visualize_iterations_total(sampler_penalty; save=true)
visualize_optimality_merit(sampler_penalty; save=true)


include("../../src/utils/monte_carlo_visualization_latex.jl")
axis_H_cond = visualize_latex_H_cond(sampler_direct; save=true)
axis_solve_time = visualize_latex_solve_time(sampler_direct; save=true)
axis_cmax = visualize_latex_cmax(sampler_direct; save=true)
axis_iterations_total = visualize_latex_iterations_total(sampler_direct; save=true)
axis_optimality_merit = visualize_latex_optimality_merit(sampler_direct; save=true)

include("../../src/utils/monte_carlo_visualization_latex.jl")
axis_H_cond = visualize_latex_H_cond(sampler_penalty; save=true)
axis_solve_time = visualize_latex_solve_time(sampler_penalty; save=true)
axis_cmax = visualize_latex_cmax(sampler_penalty; save=true)
axis_iterations_total = visualize_latex_iterations_total(sampler_penalty; save=true)
axis_optimality_merit = visualize_latex_optimality_merit(sampler_penalty; save=true)

gp = visualize_latex_sampler(sampler_direct; save=true)
gp = visualize_latex_sampler(sampler_penalty; save=true)



a = 100
a = 100
a = 100
a = 100
a = 100


cmax_direct = [[sampler_direct.stats.cmax[i], i] for i in 1:sampler_direct.stats.iterations]
sort(cmax_direct)
cmax_viol_direct = []
for i = 1:sampler_direct.stats.iterations
    if sampler_direct.stats.cmax[i] > 1e-3
        push!(cmax_viol_direct, i)
    end
end
cmax_viol_direct

opt_direct = [[sampler_direct.stats.optimality_merit[i], i] for i in 1:sampler_direct.stats.iterations]
sort(opt_direct)
opt_viol_direct = []
for i = 1:sampler_direct.stats.iterations
    if sampler_direct.stats.optimality_merit[i] > 1e-2
        push!(opt_viol_direct, i)
    end
end
opt_viol_direct








cmax_penalty = [[sampler_penalty.stats.cmax[i], i] for i in 1:sampler_penalty.stats.iterations]
sort(cmax_penalty)
cmax_viol_penalty = []
for i = 1:sampler_penalty.stats.iterations
    if sampler_penalty.stats.cmax[i] > 1e-3
        push!(cmax_viol_penalty, i)
    end
end
cmax_viol_penalty

opt_penalty = [[sampler_penalty.stats.optimality_merit[i], i] for i in 1:sampler_penalty.stats.iterations]
sort(opt_penalty)
opt_viol_penalty = []
for i = 1:sampler_penalty.stats.iterations
    if sampler_penalty.stats.optimality_merit[i] > 1e-2
        push!(opt_viol_penalty, i)
    end
end
opt_viol_penalty







# # using Plots
# using PGFPlots
# # using ImageMagick
# # pgfplots()
# histogram(randn(1000), bins=:scott, weights=repeat(1:5, outer=200))
#
#
#
# save("myfile.tex", a)
# PGFPlots.save("const" * ".tikz", include_preamble=false, a)
#
#
# push!(group_plot, axis)


sort(sampler_direct.stats.cmax)
sort(sampler_direct.stats.optimality_merit)[1:850]
sort(sampler_direct.stats.solve_time)[1:860]
sort(sampler_direct.stats.iterations_total)[end-14:end]
sort(sampler_penalty.stats.iterations_total)[end-14:end]


mean(sampler_direct.stats.solve_time)
median(sampler_direct.stats.solve_time)
var(sampler_direct.stats.solve_time)

mean(sampler_direct.stats.iterations_total)
median(sampler_direct.stats.iterations_total)
var(sampler_direct.stats.iterations_total)

sort(sampler_penalty.stats.cmax)

mean(sampler_penalty.stats.solve_time)
median(sampler_penalty.stats.solve_time)
var(sampler_penalty.stats.solve_time)

mean(sampler_penalty.stats.iterations_total)
median(sampler_penalty.stats.iterations_total)
var(sampler_penalty.stats.iterations_total)
