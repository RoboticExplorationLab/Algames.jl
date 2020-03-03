using ALGAMES
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
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
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [ # Initial state
               -0.50, -0.15,  0.00,  0.60, #lane1
               1.40,   0.15,  pi,    0.60, #lane2
               0.43,  -0.30,  pi/2,  0.12, #lane4
                ]
xf = @SVector [ # Final State
               1.30,  -0.15,  0.00, 0.60, #lane1
               -0.30,  0.15,  pi,   0.60, #lane2
               0.43,   0.35,  pi/2, 0.30, #lane4
               ]
dxf = @SVector [ # Variation of the fnal state when executing an MPC
               0.60,  0.00,  0.00,  0.00, #lane1
              -0.60,  0.00,  0.00,  0.00, #lane2
               0.00,  0.12,  0.00,  0.00, #lane2
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
    0., 1., 1., 1.]
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
actors_types = [:car, :car, :pedestrian] # For visualization
top_road_length = 4.0
road_width = 0.60
bottom_road_length = 1.0
cross_width = 0.25
bound_radius = 0.05
lanes = [1, 2, 5] # Lane indices for each player
scenario = TIntersectionScenario(
    top_road_length, road_width, bottom_road_length,
    cross_width, car_radii, actors_types, bound_radius)

# Create constraints
conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(conSet, car_radii, px, p, con_inds)
# Add scenario specific constraints
add_scenario_constraints(conSet, scenario, lanes, px, con_inds; constraint_type=:constraint)


prob = GameProblem(model, obj, conSet, x0, xf, Z, N, tf)
opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
	log_level=Logging.Debug)
solver = DirectGamesSolver(prob, opts)
reset!(solver, reset_type=:full)
@time solve!(solver)

solver.opts.log_level = AG.Logging.Warn
@btime timing_solve(solver)

# reset!(solver, reset_type=:full)
# @profiler AG.solve!(solver)

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

vis=Visualizer()
anim=MeshCat.Animation()
open(vis)
# Execute this line after the MeshCat tab is open
vis, anim = animation(solver, scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=false)
