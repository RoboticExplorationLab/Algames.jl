################################################################################
# Drone Example
################################################################################
# using Algames
using StaticArrays
using LinearAlgebra


vis = Visualizer()
open(vis)

T = Float64

# Define the dynamics of the system
p = 2 # Number of players
d = 3 # Number of dimensions
model = DoubleIntegratorGame(p=p, d=d) # game with 3 players with double integrator dynamics in 3D
n = model.n
m = model.m

# Define the horizon of the problem
N = 20 # N time steps
dt = 0.1 # each step lasts 0.1 second
probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

# Define the objective of each player
# We use a LQR cost
Q = [Diagonal(10*ones(SVector{model.ni[i],T})) for i=1:p] # Quadratic state cost
R = [Diagonal(0.1*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
# Desrired state
xf = [SVector{model.ni[1],T}([2,+0.7,+0.3,0,0,0]),
      SVector{model.ni[2],T}([2,-0.7,-0.3,0,0,0]),
      ]
# Desired control
uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
# Objectives of the game
game_obj = GameObjective(Q,R,xf,uf,N,model)
radius = 1.0*ones(p)
μ = 5.0*ones(p)
add_collision_cost!(game_obj, radius, μ)

# Define the constraints that each player must respect
game_con = GameConstraintValues(probsize)
# Add collision avoidance
radius = 0.08
add_spherical_collision_avoidance!(game_con, radius)
# Add control bounds
u_max =  5*ones(SVector{m,T})
u_min = -5*ones(SVector{m,T})
add_control_bound!(game_con, u_max, u_min)
# # Add state bounds for player 1
x_max =  5*ones(SVector{n,T})
x_min = -5*ones(SVector{n,T})
add_state_bound!(game_con, 1, x_max, x_min)
# Add wall constraint
# walls = [Wall3D([1.,0,0], [2.,0,0], [2.,1,0], [0.,0,1])]
room_walls = [
    # Wall3D([1.,0.5,0], [2.,0,1], [2.5,1.0,1], [0.,0,1]),
    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00,  1.00, -1.00], [0.00,  0.00, -1.00]),
    Wall3D([-3.00, -1.00,  1.00], [3.00, -1.00,  1.00], [3.00,  1.00,  1.00], [0.00,  0.00,  1.00]),
    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00, -1.00,  1.00], [0.00, -1.00,  0.00]),
    Wall3D([-3.00,  1.00, -1.00], [3.00,  1.00, -1.00], [3.00,  1.00,  1.00], [0.00,  1.00,  0.00]),
    ]
door_walls = [
    Wall3D([-0.75, -1.00,  1.00], [0.00, -1.00,  0.25], [0.00,  1.00,  0.25], [ sqrt(2),  0.00,  sqrt(2)]),
    Wall3D([ 0.75, -1.00,  1.00], [0.00, -1.00,  0.25], [0.00,  1.00,  0.25], [-sqrt(2),  0.00,  sqrt(2)]),

    Wall3D([-0.75, -1.00, -1.00], [0.00, -1.00, -0.25], [0.00,  1.00, -0.25], [ sqrt(2),  0.00, -sqrt(2)]),
    Wall3D([ 0.75, -1.00, -1.00], [0.00, -1.00, -0.25], [0.00,  1.00, -0.25], [-sqrt(2),  0.00, -sqrt(2)]),
    #
    Wall3D([-0.75, -1.00, -1.00], [0.00, -0.25, -1.00], [0.00, -0.25,  1.00], [ sqrt(2), -sqrt(2),  0.00]),
    Wall3D([ 0.75, -1.00, -1.00], [0.00, -0.25, -1.00], [0.00, -0.25,  1.00], [-sqrt(2), -sqrt(2),  0.00]),
    #
    Wall3D([-0.75,  1.00, -1.00], [0.00,  0.25, -1.00], [0.00,  0.25,  1.00], [ sqrt(2),  sqrt(2),  0.00]),
    Wall3D([ 0.75,  1.00, -1.00], [0.00,  0.25, -1.00], [0.00,  0.25,  1.00], [-sqrt(2),  sqrt(2),  0.00]),
    ]
add_wall_constraint!(game_con, room_walls)
add_wall_constraint!(game_con, door_walls)

build_wall!(vis, room_walls, α=0.1, name=:Room)
build_wall!(vis, door_walls, α=0.1, name=:Door1)






# Define the initial state of the system
x0 = SVector{model.n,T}([
   -1.2,-1.0,
   -0.4, 0.4,
    0.6, 0.6,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    ])

# Define the Options of the solver
opts = Options()
# Define the game problem
prob = GameProblem(N,dt,x0,model,opts,game_obj,game_con)

# Solve the problem
@time newton_solve!(prob)

# Visualize the Results
const Algames = Main
using Plots
Algames.plot_traj!(prob.model, prob.pdtraj.pr)
Algames.plot_violation!(prob.stats)

build_traj!(vis, model, prob.pdtraj.pr, α=1.0, name=:Traj)

# build_robot!(vis, model, name=:Robot, α=1.0)
# set_robot!(vis, model, TrajectoryOptimization.state(prob.pdtraj.pr[1]), name=:Robot)
# animate_robot!(vis, model, α=1.4, name=:Robot)
anim = visualize_robot!(vis, model, prob.pdtraj.pr)
