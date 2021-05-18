################################################################################
# Drone Example
################################################################################
# using Algames
using StaticArrays
using LinearAlgebra
using MeshCat

vis  = Visualizer()
open(vis)
T = Float64

# Define the dynamics of the system
p = 4 # Number of players
model  = QuadrotorGame(p=p) # game with 3 players with double integrator dynamics in 3D
n = model.n
m = model.m

# Define the horizon of the problem
N = 20 # N time steps
dt = 0.40 # each step lasts 0.1 second
probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

# Define the objective of each player
# We use a LQR cost
Q = [Diagonal(5*SVector{model.ni[i],T}([1*[0,0.3,0.3,1,1,1]; 5*[1,0.3,0.3,1,1,1]])) for i=1:p] # Quadratic state cost
R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
# Desrired state
xf = [SVector{model.ni[1],T}([1.0,-0.5,+0.3,0,0,0, 0.3,0,0, 0,0,0]),
	  SVector{model.ni[2],T}([1.5,-0.7,-0.3,0,0,0, 0.3,0,0, 0,0,0]),
	  SVector{model.ni[3],T}([0.8,-0.0,+0.5,0,0,0, 0.3,0,0, 0,0,0]),
      SVector{model.ni[3],T}([0.8,-0.0,-0.5,0,0,0, 0.3,0,0, 0,0,0]),
      ]
xf = xf[1:p]
# Desired control
uf = [- model.mass * model.gravity[end] / 4 / model.kf * ones(SVector{model.mi[i],T}) for i=1:p]
# Objectives of the game
game_obj = GameObjective(Q,R,xf,uf,N,model)
radius = 0.2*ones(p)
μ = 20.0*ones(p)
add_collision_cost!(game_obj, radius, μ)

# Define the constraints that each player must respect
game_con = GameConstraintValues(probsize)
# Add collision avoidance
radius = 0.08
add_spherical_collision_avoidance!(game_con, radius)
# Add wall constraint
room_walls = [
    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00,  1.00, -1.00], [0.00,  0.00, -1.00]),
    Wall3D([-3.00, -1.00,  1.00], [3.00, -1.00,  1.00], [3.00,  1.00,  1.00], [0.00,  0.00,  1.00]),
    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00, -1.00,  1.00], [0.00, -1.00,  0.00]),
    Wall3D([-3.00,  1.00, -1.00], [3.00,  1.00, -1.00], [3.00,  1.00,  1.00], [0.00,  1.00,  0.00]),
    ]
door_cylinders = [
	CylinderWall([0.00, -1.00, -1.00], :z, 2.00, 0.95),
	CylinderWall([0.00,  1.00, -1.00], :z, 2.00, 0.95),
	CylinderWall([0.00, -1.00, -1.00], :y, 2.00, 0.95),
	CylinderWall([0.00, -1.00,  1.00], :y, 2.00, 0.95),
	]
add_wall_constraint!(game_con, room_walls)
add_wall_constraint!(game_con, door_cylinders)

build_wall!(vis, room_walls, α=0.1, name=:Room)
build_cylinder!(vis, door_cylinders, radius=radius, α=0.1, name=:Door1)
a = ones(4)
b = zeros(4)

reshape(vcat(a', b'), (8))

reshape([a'; b'], (8,1))
# Define the initial state of the system
x0 = [
	[-1.0, -0.4, -0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	[-0.8,  0.4,  0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	[-1.4,  0.2,  0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	[-1.6, -0.3, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	]
x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))

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
build_xf!(vis, model, xf, α=1.0, name=:Xf)

anim = visualize_robot!(vis, model, prob.pdtraj.pr)

plot(hcat([state(prob.pdtraj.pr[i]) for i=1:N]...)')
plot(hcat([control(prob.pdtraj.pr[i]) for i=1:N]...)')
