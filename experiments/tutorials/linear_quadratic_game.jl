using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Instantiate dynamics model
model = DoubleIntegratorGame(p=2)
n,m,pu,p = size(model)
T = Float64

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [-0.50, -0.50, # x
			    0.10, -0.10, # y
                0.50,  0.40, # xdot
				0.00,  0.00] # ydot
xf = @SVector [ 0.50,  0.50, # x
			   -0.10,  0.10, # y
                0.40,  0.30, # xdot
				0.00,  0.10] # ydot

# Define a quadratic cost for each player
diag_Q = [SVector{n}([ 1.0, -0.1,  1.0, -0.1,  1.0, -0.1,  1.0, -0.1]), 	# Player 1 cost
	      SVector{n}([-0.1,  1.0, -0.1,  1.0, -0.1,  1.0, -0.1,  1.0])] 	# Player 2 cost
Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[i]))) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N,checks=false) for i=1:p]

# Build problem
actor_radius = 0.08
actors_radii = [actor_radius for i=1:p]

# Create constraints
conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Define the problem
prob = GameProblem(model, obj, xf, tf, constraints=conSet, x0=x0, N=N)
opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
	log_level=AG.Logging.Warn)
solver = DirectGamesSolver(prob, opts)
@time solve!(solver)

# X = TO.states(solver)
# U = TO.controls(solver)
# visualize_state(X)
# visualize_control(U,pu)
# visualize_trajectory_car(solver)
# visualize_collision_avoidance(solver)
# visualize_dynamics(solver)
# visualize_optimality_merit(solver)
# visualize_α(solver)
# visualize_cmax(solver)
