using ALGAMES
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
using StatsBase
const TO = TrajectoryOptimization
const AG = ALGAMES

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
N = 41   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [
               -0.50, -0.15,  0.00, 0.60, #lane1
                1.40,  0.15,  pi,   0.60, #lane2
                0.43, -0.30,  pi/2, 0.20, #lane4
                ]
xf = @SVector [
                1.30, -0.15,  0.00, 0.60, #lane1
               -0.30,  0.15,  pi,   0.60, #lane2
                0.43,  0.35,  pi/2, 0.20, #lane4
               ]

# Define the movement of the goal state xf with time
# xf <- xf + dxf * dt
# This allows the goal of each vehicle to be updated as they
# are solving the MPC. We move the goal forward to keep
# to incentivize the vehicles to continue moving forward.
dxf = @SVector [
               0.60,  0.00,  0.00,  0.00, #lane1
              -0.60,  0.00,  0.00,  0.00, #lane2
               0.00,  0.20,  0.00,  0.00, #lane4
              ]

# Define a quadratic cost
diag_Q1 = @SVector [ # Player 1 state cost
	0., 1., 1., 1.,
	0., 0., 0., 0.,
	0., 0., 0., 0.]
diag_Q2 = @SVector [ # Player 2 state cost
	0., 0., 0., 0.,
	0., 1., 1., 1.,
	0., 0., 0., 0.]
diag_Q3 = @SVector [ # Player 3 state cost
	0., 0., 0., 0.,
	0., 0., 0., 0.,
	1., 0., 1., 1.]
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
pr = 3
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
actors_radii_prog = [[actor_radius + (j-1)*0.01 for i=1:p] for j=1:pr]

con_inds_prog = [j:j for j=1:pr]
con_inds_prog[end] = pr:N

actors_types = [:car, :car, :pedestrian]
top_road_length = 4.0
road_width = 0.60 #0.60
bottom_road_length = 1.0
cross_width = 0.25
bound_radius = 0.05
lanes = [1, 2, 5]
scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length, cross_width,
    actors_radii, actors_types, bound_radius)

# Create constraint sets
algames_conSet = ConstraintSet(n,m,N)
# Add collision avoidance constraints
for j = 1:pr
	add_collision_avoidance(algames_conSet, actors_radii_prog[j], px, p, con_inds_prog[j])
end
# Add scenario specific constraints (road boundaries)
con_inds = 2:N # Indices where the constraints will be applied
add_scenario_constraints(algames_conSet, scenario, lanes, px, con_inds; constraint_type=:constraint)

algames_prob = GameProblem(model, obj, algames_conSet, x0, xf, Z, N, tf)
algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
	min_steps_per_iteration=0,
	record_condition=false,
    inner_iterations=20,
    iterations_linesearch=10,
    log_level=ALGAMES.Logging.Warn)
algames_solver = DirectGamesSolver(algames_prob, algames_opts)
@time solve!(algames_solver)


state_noise = 5. * SVector{n}([
	0.008, 0.008, 2*pi/72, 0.03, #+-50cm, +-50cm, +-25deg, +-12.5% per second
	0.008, 0.008, 2*pi/72, 0.03,
	0.008, 0.008, 2*pi/72, 0.03])
opts_mpc = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=1000,
	N_mpc=50,
	mpc_tf=4.0,
	min_δt=0.005,
    dxf=dxf,
	noise=state_noise)
mpc_solver = MPCGamesSolver(algames_solver, opts_mpc)
reset!(mpc_solver, reset_type=:full)
solve!(mpc_solver; wait=false)
resample!(mpc_solver)


vis=AG.Visualizer()
anim=AG.MeshCat.Animation()
open(vis)
anim_opts = AnimationOptions(display_actors=true,
	display_trajectory=true)
# Execute this line after the MeshCat tab is open
vis, anim = animation(mpc_solver, scenario;
	vis=vis, anim=anim, opts=anim_opts)

iter = mpc_solver.stats.iterations
mean_solve_time = mean(mpc_solver.stats.solve_time[1:iter])
update_freq = mpc_solver.stats.iterations/mpc_solver.stats.time
std_solve_time = StatsBase.std(mpc_solver.stats.solve_time)
largest_solve_time = maximum(mpc_solver.stats.solve_time)


# To evaluate the update frequency of the MPC solver we average over many samples.
samples = 100
times = zeros(0)
cmax = zeros(samples)
mpc_solver.opts.log_level = AG.Logging.Debug
for k = 1:samples
	@show k
    algames_solver = DirectGamesSolver(algames_prob, algames_opts)
    mpc_solver = MPCGamesSolver(algames_solver, opts_mpc)

    reset!(mpc_solver, reset_type=:full)
    solve!(mpc_solver, wait=false)
    i = mpc_solver.stats.iterations
    cmax[k] = maximum(mpc_solver.stats.cmax[1:i])
    push!(times, mpc_solver.stats.solve_time[1:i]...)
end

# Average MPC frequency
freq = length(times) / sum(times) # 176 Hz
# Mean solve time
mean_solve_time = sum(times) / length(times) #
# Maximum constraint violation across samples.
# This constraint violation should be comparable to the tolerance
# we used in the solver.
max_constraint_violation = maximum(cmax)
solver_violation_tolerance = mpc_solver.solver.opts.constraint_tolerance
# If the max_constraint_violation <= solver_violation_tolerance
# this means that all the open-loop plans have converged to constraint
# satisfaction during the MPC solve.
