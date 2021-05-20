################################################################################
# Drone Example
################################################################################
# using Algames
using StaticArrays
using LinearAlgebra
using MeshCat


function quadrotor_experiment(vis::Visualizer; p::Int=2, animate::Bool=false, S::Int=1)
	T = Float64
	# Define the dynamics of the system
	model  = QuadrotorGame(p=p) # game with 3 players with double integrator dynamics in 3D

	# Define the horizon of the problem
	N = 20 # N time steps
	dt = 0.40 # each step lasts 0.1 second
	probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

	# Define the objective of each player
	# We use a LQR cost
	Q = [Diagonal(5*SVector{model.ni[i],T}([1*[0,0.3,0.3,1,1,1]; 5*[1,0.3,0.3,1,1,1]])) for i=1:p] # Quadratic state cost
	R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
	# Desrired state
	xf = [SVector{12,T}([1.0,-0.5,+0.3,0,0,0, 0.3,0,0, 0,0,0]),
		  SVector{12,T}([1.5,-0.7,-0.3,0,0,0, 0.3,0,0, 0,0,0]),
		  SVector{12,T}([0.8, 0.1,+0.5,0,0,0, 0.3,0,0, 0,0,0]),
	      SVector{12,T}([0.8,-0.1,-0.5,0,0,0, 0.3,0,0, 0,0,0]),
	      ]
	xf = xf[1:p]
	# Desired control
	uf = [- model.mass * model.gravity[end] / 4 / model.kf * ones(SVector{model.mi[i],T}) for i=1:p]
	# Objectives of the game
	game_obj = GameObjective(Q,R,xf,uf,N,model)
	radius = 0.5*ones(p)
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
	t_elapsed = @elapsed for i = 1:S
		newton_solve!(prob)
	end
	t_elapsed /= S

	# Visualize the Results
	plot_traj!(prob.model, prob.pdtraj.pr)
	plot_violation!(prob.stats)

	build_traj!(vis, model, prob.pdtraj.pr, α=1.0, name=:Traj)
	build_xf!(vis, model, xf, α=1.0, name=:Xf)

	animate && (anim = visualize_robot!(vis, model, prob.pdtraj.pr))
	return prob, t_elapsed
end

function double_integrator_3D_experiment(vis::Visualizer; p::Int=2, animate::Bool=false, S::Int=1)
	T = Float64
	# Define the dynamics of the system
	model  = DoubleIntegratorGame(p=p, d=3) # game with 3 players with double integrator dynamics in 3D

	# Define the horizon of the problem
	N = 20 # N time steps
	dt = 0.40 # each step lasts 0.1 second
	probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

	# Define the objective of each player
	# We use a LQR cost
	Q = [Diagonal(5*SVector{model.ni[i],T}([1*[0,0.3,0.3]; 5*[1,0.3,0.3]])) for i=1:p] # Quadratic state cost
	R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
	# Desrired state
	xf = [SVector{6,T}([1.0,-0.5,+0.3, 0.3,0,0]),
		  SVector{6,T}([1.5,-0.7,-0.3, 0.3,0,0]),
		  SVector{6,T}([0.8, 0.1,+0.5, 0.3,0,0]),
	      SVector{6,T}([0.8,-0.1,-0.5, 0.3,0,0]),
	      ]
	xf = xf[1:p]
	# Desired control
	uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
	# Objectives of the game
	game_obj = GameObjective(Q,R,xf,uf,N,model)
	radius = 0.5*ones(p)
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

	# Define the initial state of the system
	x0 = [
		[-1.0, -0.4, -0.6, 0.0, 0.0, 0.0],
		[-0.8,  0.4,  0.6, 0.0, 0.0, 0.0],
		[-1.4,  0.2,  0.3, 0.0, 0.0, 0.0],
		[-1.6, -0.3, -0.3, 0.0, 0.0, 0.0],
		]
	x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))

	# Define the Options of the solver
	opts = Options()
	# Define the game problem
	prob = GameProblem(N,dt,x0,model,opts,game_obj,game_con)

	# Solve the problem
	t_elapsed = @elapsed for i = 1:S
		newton_solve!(prob)
	end
	t_elapsed /= S

	# Visualize the Results
	plot_traj!(prob.model, prob.pdtraj.pr)
	plot_violation!(prob.stats)

	build_traj!(vis, model, prob.pdtraj.pr, α=1.0, name=:Traj)
	build_xf!(vis, model, xf, α=1.0, name=:Xf)

	animate && (anim = visualize_robot!(vis, model, prob.pdtraj.pr))
	return prob, t_elapsed
end

function double_integrator_2D_experiment(vis::Visualizer; p::Int=2, animate::Bool=false, S::Int=1)
	T = Float64
	# Define the dynamics of the system
	model  = DoubleIntegratorGame(p=p, d=2) # game with 3 players with double integrator dynamics in 3D

	# Define the horizon of the problem
	N = 20 # N time steps
	dt = 0.40 # each step lasts 0.1 second
	probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

	# Define the objective of each player
	# We use a LQR cost
	Q = [Diagonal(5*SVector{model.ni[i],T}([1*[0,0.3]; 5*[1,0.3]])) for i=1:p] # Quadratic state cost
	R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
	# Desrired state
	xf = [SVector{4,T}([1.0,-0.5, 0.3,0]),
		  SVector{4,T}([1.5,-0.7, 0.3,0]),
		  SVector{4,T}([0.8, 0.1, 0.3,0]),
	      SVector{4,T}([0.8,-0.1, 0.3,0]),
	      ]
	xf = xf[1:p]
	# Desired control
	uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
	# Objectives of the game
	game_obj = GameObjective(Q,R,xf,uf,N,model)
	radius = 0.5*ones(p)
	μ = 20.0*ones(p)
	add_collision_cost!(game_obj, radius, μ)

	# Define the constraints that each player must respect
	game_con = GameConstraintValues(probsize)
	# Add collision avoidance
	radius = 0.08
	add_spherical_collision_avoidance!(game_con, radius)
	# Add wall constraint
	room_walls = [
	    Wall([-3.00, -1.00], [3.00, -1.00], [0.00, -1.00]),
	    Wall([-3.00,  1.00], [3.00,  1.00], [0.00,  1.00]),
	    ]
	xc = [ 0.00, 0.00]
	yc = [-1.00, 1.00]
	circle_radius = [0.95, 0.95]

	add_wall_constraint!(game_con, room_walls)
	add_circle_constraint!(game_con, xc, yc, circle_radius)

	room_walls_vis = [
	    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00, -1.00,  1.00], [0.00, -1.00,  0.00]),
	    Wall3D([-3.00,  1.00, -1.00], [3.00,  1.00, -1.00], [3.00,  1.00,  1.00], [0.00,  1.00,  0.00]),
	    ]
	door_cylinders_vis = [
		CylinderWall([0.00, -1.00, -1.00], :z, 2.00, 0.95),
		CylinderWall([0.00,  1.00, -1.00], :z, 2.00, 0.95),
		]
	build_wall!(vis, room_walls_vis, α=0.1, name=:Room)
	build_cylinder!(vis, door_cylinders_vis, radius=radius, α=0.1, name=:Door1)

	# Define the initial state of the system
	x0 = [
		[-1.0, -0.4, 0.0, 0.0],
		[-0.8,  0.4, 0.0, 0.0],
		[-1.4,  0.2, 0.0, 0.0],
		[-1.6, -0.3, 0.0, 0.0],
		]
	x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))

	# Define the Options of the solver
	opts = Options()
	# Define the game problem
	prob = GameProblem(N,dt,x0,model,opts,game_obj,game_con)

	# Solve the problem
	t_elapsed = @elapsed for i = 1:S
		newton_solve!(prob)
	end
	t_elapsed /= S

	# Visualize the Results
	plot_traj!(prob.model, prob.pdtraj.pr)
	plot_violation!(prob.stats)

	build_traj!(vis, model, prob.pdtraj.pr, α=1.0, name=:Traj)
	build_xf!(vis, model, xf, α=1.0, name=:Xf)

	animate && (anim = visualize_robot!(vis, model, prob.pdtraj.pr))
	return prob, t_elapsed
end

function unicycle_experiment(vis::Visualizer; p::Int=2, animate::Bool=false, S::Int=1)
	T = Float64
	# Define the dynamics of the system
	model = UnicycleGame(p=p) # game with 3 players with double integrator dynamics in 3D

	# Define the horizon of the problem
	N = 20 # N time steps
	dt = 0.40 # each step lasts 0.1 second
	probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

	# Define the objective of each player
	# We use a LQR cost
	Q = [Diagonal(5*SVector{model.ni[i],T}([1*[0,0.3]; 5*[0.1,1.0]])) for i=1:p] # Quadratic state cost
	R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
	# Desrired state
	xf = [SVector{4,T}([1.0,-0.5, 0.0, 0.3]),
		  SVector{4,T}([1.5,-0.7, 0.0, 0.3]),
		  SVector{4,T}([0.8, 0.1, 0.0, 0.3]),
	      SVector{4,T}([0.8,-0.1, 0.0, 0.3]),
	      ]
	xf = xf[1:p]
	# Desired control
	uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
	# Objectives of the game
	game_obj = GameObjective(Q,R,xf,uf,N,model)
	radius = 0.5*ones(p)
	μ = 20.0*ones(p)
	add_collision_cost!(game_obj, radius, μ)

	# Define the constraints that each player must respect
	game_con = GameConstraintValues(probsize)
	# Add collision avoidance
	radius = 0.08
	add_spherical_collision_avoidance!(game_con, radius)
	# Add wall constraint
	room_walls = [
	    Wall([-3.00, -1.00], [3.00, -1.00], [0.00, -1.00]),
	    Wall([-3.00,  1.00], [3.00,  1.00], [0.00,  1.00]),
	    ]
	xc = [ 0.00, 0.00]
	yc = [-1.00, 1.00]
	circle_radius = [0.95, 0.95]

	add_wall_constraint!(game_con, room_walls)
	add_circle_constraint!(game_con, xc, yc, circle_radius)

	room_walls_vis = [
	    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00, -1.00,  1.00], [0.00, -1.00,  0.00]),
	    Wall3D([-3.00,  1.00, -1.00], [3.00,  1.00, -1.00], [3.00,  1.00,  1.00], [0.00,  1.00,  0.00]),
	    ]
	door_cylinders_vis = [
		CylinderWall([0.00, -1.00, -1.00], :z, 2.00, 0.95),
		CylinderWall([0.00,  1.00, -1.00], :z, 2.00, 0.95),
		]
	build_wall!(vis, room_walls_vis, α=0.1, name=:Room)
	build_cylinder!(vis, door_cylinders_vis, radius=radius, α=0.1, name=:Door1)

	# Define the initial state of the system
	x0 = [
		[-1.0, -0.4, 0.0, 0.0],
		[-0.8,  0.4, 0.0, 0.0],
		[-1.4,  0.2, 0.0, 0.0],
		[-1.6, -0.3, 0.0, 0.0],
		]
	x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))

	# Define the Options of the solver
	opts = Options()
	# Define the game problem
	prob = GameProblem(N,dt,x0,model,opts,game_obj,game_con)

	# Solve the problem
	t_elapsed = @elapsed for i = 1:S
		newton_solve!(prob)
	end
	t_elapsed /= S

	# Visualize the Results
	plot_traj!(prob.model, prob.pdtraj.pr)
	plot_violation!(prob.stats)

	build_traj!(vis, model, prob.pdtraj.pr, α=1.0, name=:Traj)
	build_xf!(vis, model, xf, α=1.0, name=:Xf)

	animate && (anim = visualize_robot!(vis, model, prob.pdtraj.pr))
	return prob, t_elapsed
end

function bicycle_experiment(vis::Visualizer; p::Int=2, animate::Bool=false, S::Int=1)
	T = Float64
	# Define the dynamics of the system
	model = BicycleGame(p=p) # game with 3 players with double integrator dynamics in 3D

	# Define the horizon of the problem
	N = 20 # N time steps
	dt = 0.40 # each step lasts 0.1 second
	probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

	# Define the objective of each player
	# We use a LQR cost
	Q = [Diagonal(5*SVector{model.ni[i],T}([1*[0,0.3]; 5*[1.0, 0.1]])) for i=1:p] # Quadratic state cost
	R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
	# Desrired state
	xf = [SVector{4,T}([1.0,-0.5, 0.3, 0.0]),
		  SVector{4,T}([1.5,-0.7, 0.3, 0.0]),
		  SVector{4,T}([0.8, 0.1, 0.3, 0.0]),
	      SVector{4,T}([0.8,-0.1, 0.3, 0.0]),
	      ]
	xf = xf[1:p]
	# Desired control
	uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
	# Objectives of the game
	game_obj = GameObjective(Q,R,xf,uf,N,model)
	radius = 0.5*ones(p)
	μ = 20.0*ones(p)
	add_collision_cost!(game_obj, radius, μ)

	# Define the constraints that each player must respect
	game_con = GameConstraintValues(probsize)
	# Add collision avoidance
	radius = 0.08
	add_spherical_collision_avoidance!(game_con, radius)
	# Add wall constraint
	room_walls = [
	    Wall([-3.00, -1.00], [3.00, -1.00], [0.00, -1.00]),
	    Wall([-3.00,  1.00], [3.00,  1.00], [0.00,  1.00]),
	    ]
	xc = [ 0.00, 0.00]
	yc = [-1.00, 1.00]
	circle_radius = [0.95, 0.95]

	add_wall_constraint!(game_con, room_walls)
	add_circle_constraint!(game_con, xc, yc, circle_radius)

	room_walls_vis = [
	    Wall3D([-3.00, -1.00, -1.00], [3.00, -1.00, -1.00], [3.00, -1.00,  1.00], [0.00, -1.00,  0.00]),
	    Wall3D([-3.00,  1.00, -1.00], [3.00,  1.00, -1.00], [3.00,  1.00,  1.00], [0.00,  1.00,  0.00]),
	    ]
	door_cylinders_vis = [
		CylinderWall([0.00, -1.00, -1.00], :z, 2.00, 0.95),
		CylinderWall([0.00,  1.00, -1.00], :z, 2.00, 0.95),
		]
	build_wall!(vis, room_walls_vis, α=0.1, name=:Room)
	build_cylinder!(vis, door_cylinders_vis, radius=radius, α=0.1, name=:Door1)

	# Define the initial state of the system
	x0 = [
		[-1.0, -0.4, 0.0, 0.0],
		[-0.8,  0.4, 0.0, 0.0],
		[-1.4,  0.2, 0.0, 0.0],
		[-1.6, -0.3, 0.0, 0.0],
		]
	x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))

	# Define the Options of the solver
	opts = Options()
	# Define the game problem
	prob = GameProblem(N,dt,x0,model,opts,game_obj,game_con)

	# Solve the problem
	t_elapsed = @elapsed for i = 1:S
		newton_solve!(prob)
	end
	t_elapsed /= S

	# Visualize the Results
	plot_traj!(prob.model, prob.pdtraj.pr)
	plot_violation!(prob.stats)

	build_traj!(vis, model, prob.pdtraj.pr, α=1.0, name=:Traj)
	build_xf!(vis, model, xf, α=1.0, name=:Xf)

	animate && (anim = visualize_robot!(vis, model, prob.pdtraj.pr))
	return prob, t_elapsed
end

S = 100

# vis1  = Visualizer()
# open(vis1)
prob1 = []
t1 = []
for p = 1:4
	prob_i, t_i = quadrotor_experiment(vis1, p=p, animate=false, S=S)
	push!(prob1, prob_i)
	push!(t1, t_i)
end
@show t1
# t = [0.12, 0.45, 1.34, 3.32]

# vis2  = Visualizer()
# open(vis2)
prob2 = []
t2 = []
for p = 1:4
	prob_i, t_i = double_integrator_3D_experiment(vis2, p=p, animate=false, S=S)
	push!(prob2, prob_i)
	push!(t2, t_i)
end
@show t2
# t = [0.0159, 0.214, 0.9569, 0.9119]

# vis3  = Visualizer()
# open(vis3)
prob3 = []
t3 = []
for p = 1:4
	prob_i, t_i = bicycle_experiment(vis3, p=p, animate=false, S=S)
	push!(prob3, prob_i)
	push!(t3, t_i)
end
@show t3
# t = [0.029, 0.125, 0.26, 0.454]

# vis4  = Visualizer()
# open(vis4)
prob4 = []
t4 = []
for p = 1:4
	prob_i, t_i = unicycle_experiment(vis4, p=p, animate=false, S=S)
	push!(prob4, prob_i)
	push!(t4, t_i)
end
@show t4
# t = [0.029, 0.132, 0.276, 0.745]

# vis5  = Visualizer()
# open(vis5)
prob5 = []
t5 = []
for p = 1:4
	prob_i, t_i = double_integrator_2D_experiment(vis5, p=p, animate=false, S=S)
	push!(prob5, prob_i)
	push!(t5, t_i)
end
@show t5
# t = [0.0095, 0.0695, 0.1640, 0.6938]

@show scn.(t1)
@show scn.(t2)
@show scn.(t3)
@show scn.(t4)
@show scn.(t5)
