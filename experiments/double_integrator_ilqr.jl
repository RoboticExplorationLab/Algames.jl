using TrajectoryOptimization
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using Plots
const TO = TrajectoryOptimization
import TrajectoryOptimization: dynamics, RK3, AbstractModel, KnotPoint, Traj, StaticBoundConstraint, GoalConstraint,
    ConstraintVals, ConstraintSets, StaticProblem, StaticALSolver
import TrajectoryOptimization: dynamics

# Define dynamics model
#   must inherit from AbstractModel
struct DoubleIntegrator{T} <: AbstractModel
    n::Int
    m::Int
    mp::T
end

DoubleIntegrator() = DoubleIntegrator(4, 2, 1.0)

#   must define this
Base.size(::DoubleIntegrator) = 4,2

#   must define this
function dynamics(model::DoubleIntegrator, x, u)
    mp = model.mp  # mass of the point mass in kg (10)
    q = x[ @SVector [1,2] ]
    qd = x[ @SVector [3,4] ]
    qdd = u/mp
    return [qd; qdd]
end

# Instantiate dynamics model
model = DoubleIntegrator()
n,m = size(model)
# n,m = size(model)

# Discretization info
tf = 5.0  # final time
N = 101   # number of knot points
dt = tf / (N-1)

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [1., 1.,0., 0.]
xf = @SVector [0., 0., 0., 0.] # (ie, swing up)

# Define a quadratic cost
Q = 1.0e-2*Diagonal(@SVector ones(n))
Qf = 100.0*Diagonal(@SVector ones(n))
R = 1.0e-1*Diagonal(@SVector ones(m))
obj = LQRObjective(Q,R,Qf,xf,N)


# Constraints
#   bound constraint
u_bnd = 0.3
bnd = StaticBoundConstraint(n,m, u_min=-u_bnd*(@SVector ones(m)), u_max=u_bnd*(@SVector ones(m)))

#   goal constraint
goal = GoalConstraint(SVector{n}(xf))

#   create constraint value types (pre-allocates arrays for values, Jacobians, active set, etc.)
#   first argument is the constraint, second is the range of knot points to which the constraint is applied
#   TODO: do this automatically when building ConstraintSets
con_bnd = ConstraintVals(bnd, 1:N-1)
con_goal = ConstraintVals(goal, N:N)

#   create constraint set for problem
# conSet = ConstraintSets([con_bnd, con_goal], N)
# conSet = ConstraintSets([con_goal], N)
conSet = ConstraintSets([], N)

# Define the initial trajectory
#  set initial states the NaNs since it will overwitten TODO: do this automatically
u0 = @SVector [0.01, 0.01]
U0 = [copy(u0) for k = 1:N]
xs = NaN*@SVector zeros(n)
us = SVector{m}(u0)
Z = [KnotPoint(xs,us,dt) for k = 1:N]
Z[end] = KnotPoint(xs,m)

# Build problem
#   TODO: better constructors
prob = StaticProblem(model, obj, conSet, x0, xf,
    deepcopy(Z), deepcopy(Z), N, dt, dt*(N-1))

# Build the solver
# opts = StaticALSolverOptions{Float64}()
# solver = StaticALSolver(prob, opts)
opts_ilqr = TrajectoryOptimization.StaticiLQRSolverOptions{Float64}()
solver_ilqr = TrajectoryOptimization.StaticiLQRSolver(prob, opts_ilqr)

# Convert to Augmented Lagrangian problem
#   TODO: convert this automatically (will incur allocations)
# prob_al = TO.convertProblem(prob, solver)

# Solve
#   reset the control before the solve so you can re-run it easily
# initial_controls!(prob_al, u0)
# solve!(prob_al, solver)
initial_controls!(prob, u0)
solve!(prob, solver_ilqr)



# Analyze results
#   look in solver.stats for stats on the solve
# println("Outer loop iterations: ", solver.stats.iterations)
println("Outer loop iterations: ", solver_ilqr.stats.iterations)
#   get final constraint violation
# println("Max violation: ", max_violation(prob_al))
#   get final cost
# println("Final cost: ", cost(prob_al))
println("Final cost: ", cost(prob))
#   extract trajectories and plot
# X = state(prob_al)
# U = control(prob_al)
X = state(prob)
U = control(prob)
plot(X, 1:2)
plot(U)
#   plot cost convergence
# plot(1:solver.stats.iterations, solver.stats.cost, yscale=:log10)
plot(1:solver_ilqr.stats.iterations, solver_ilqr.stats.cost, yscale=:log10)
#   plot constraint violation
# plot(1:solver.stats.iterations, solver.stats.c_max, yscale=:log10)

# Benchmark the result
@btime begin
    # initial_controls!(prob_al, U0)
    # solve!(prob_al, solver)
    initial_controls!(prob, U0)
    solve!(prob, solver_ilqr)
end

#   there shouldn't be any allocations
# initial_controls!(prob, U0)
# @allocated solve!(prob_al, solver)
initial_controls!(prob, U0)
@allocated solve!(prob, solver_ilqr)
