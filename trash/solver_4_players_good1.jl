using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Instantiate dynamics model
model = UnicycleGame(p=4)
n,m,pu,p = size(model)
T = Float64
px = model.px

# Discretization info
tf = 3.0  # final time
N = 21    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [# p1   # p2   # p3   # p4
			 -0.80, -1.20, -1.60, -1.00, # x
		     -0.08, -0.08, -0.08, -0.31, # y
			  0.00,  0.00,  0.00, pi/12, # θ
			  0.65,  0.65,  0.65,  0.50, # v
]
xf = @SVector [# p1   # p2   # p3   # p4
			  1.00,  0.60,  0.20,  0.80, # x
			 -0.08, -0.08, -0.08, -0.08, # y
			  0.00,  0.00,  0.00,  0.00, # θ
			  0.65,  0.65,  0.65,  0.65, # v
              ]

diag_Q = [SVector{n}([0.,  0.,  0.,  0.,
					  .1,  0.,  0.,  0.,
					  1.,  0.,  0.,  0.,
					  1.,  0.,  0.,  0.]),
	      SVector{n}([0.,  0.,  0.,  0.,
		  			  0.,  .1,  0.,  0.,
					  0.,  1.,  0.,  0.,
					  0.,  1.,  0.,  0.]),
		  SVector{n}([0.,  0.,  0.,  0.,
		  			  0.,  0.,  .1,  0.,
					  0.,  0.,  1.,  0.,
					  0.,  0.,  1.,  0.]),
		  SVector{n}([0.,  0.,  0.,  0.,
					  0.,  0.,  0.,  .1,
					  0.,  0.,  0.,  1.,
					  0.,  0.,  0.,  1.]),]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]
inflated_actors_radii = [1.6*actor_radius for i=1:p]
actors_types = [:car for i=1:p]
road_length = 34.20
road_width = 0.42
ramp_length = 17.2
ramp_angle = pi/12
ramp_merging_4_players_unicycle_penalty_scenario = MergingScenario(road_length,
	road_width, ramp_length, ramp_angle, actors_radii, actors_types)

μ_penalty = 1.0
U_lim = 1.50
