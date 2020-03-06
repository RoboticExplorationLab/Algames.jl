# ALGAMES

CODECOV_TOKEN="ccaf4089-5347-4384-a6f9-c5dd3f7c999b"

<!-- ![Build Status](https://travis-ci.org/RoboticExplorationLab/TrajectoryOptimization.jl.svg?branch=master) -->
[![codecov](https://codecov.io/gh/RoboticExplorationLab/ALGAMES.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RoboticExplorationLab/ALGAMES.jl)
<!-- [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://RoboticExplorationLab.github.io/TrajectoryOptimization.jl/dev) -->


A package for solving constrained dynamic games written in Julia. Currently, the following methods are implemented with a common interface:

[ALGAMES (Augmented Lagrangian Games Theoretic Solver)](https://rexlab.stanford.edu/papers/ALGAMES.pdf): A fast solver for constrained dynamic games that features:
  * General nonlinear cost functions
  * General nonlinear state and input constraints


This package also features:
  * Several autonomous driving environments (ramp merging, intersection crossing, etc.)
  * A Model Predictive Control (MPC) implementation of ALGAMES
  * Plotting and visualization tools

![](plots/gif/merging_4pl_2x_white.gif)

All methods utilize Julia's extensive autodifferentiation capabilities via [ForwardDiff.jl](http://www.juliadiff.org/ForwardDiff.jl/) so that the user does not need to specify derivatives of dynamics, cost, or constraint functions.

## Installation
To install ALGAMES.jl, run the following from the Julia REPL:
```julia
Pkg.add(PackageSpec(name="TrajectoryOptimization", rev="v1.3"))
Pkg.add(PackageSpec(url="https://github.com/RoboticExplorationLab/ALGAMES.jl.git"))
```

## Quick Start
To run a simple example of a 2-player 2D double-integrator:
```julia
using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization

# Define the dynamics model of the game.
struct DoubleIntegratorGame{T} <: AbstractGameModel
    n::Int  # Number of states
    m::Int  # Number of controls
    mp::T   # Mass of the point mass double integrator
	pu::Vector{Vector{Int}} # Indices of the each player's controls
	px::Vector{Vector{Int}} # Indices of the each player's x and y positions
    p::Int  # Number of players
end
DoubleIntegratorGame() = DoubleIntegratorGame(
	8,
	4,
	1.0,
	[[1,2],[3,4]],
	[[1,2],[5,6]],
	2)
Base.size(::DoubleIntegratorGame) = 8,4,[[1,2],[3,4]],2 # n,m,pu,p

# Instantiate dynamics model
model = DoubleIntegratorGame()
n,m,pu,p = size(model)
T = Float64
px = model.px
function TO.dynamics(model::DoubleIntegratorGame, x, u) # Non memory allocating dynamics
	mp = model.mp  # mass of the point mass in kg (10)
    p = model.p  # number of players
    pu = model.pu  # control vector partition for each player
    q1 = x[ @SVector [1,2] ]
    qd1 = x[ @SVector [3,4] ]
    q2 = x[ @SVector [5,6] ]
    qd2 = x[ @SVector [7,8] ]
    control1 = @SVector [u[pu_ind] for pu_ind in pu[1]]
    control2 = @SVector [u[pu_ind] for pu_ind in pu[2]]
    qdd1 = control1/mp
    qdd2 = control2/mp
    return [qd1; qdd1; qd2; qdd2]
end

# Discretization info
tf = 3.0  # final time
N = 41    # number of knot points
dt = tf / (N-1) # time step duration

# Define initial and final states (be sure to use Static Vectors!)
x0 = @SVector [-0.50,  0.10,  0.50,  0.00, #player 1
							 -0.50, -0.10,  0.40,  0.00] #player 2
xf = @SVector [ 0.50, -0.10,  0.40,  0.00, # player 1
                0.50,  0.10,  0.50,  0.80] # player 2


# Define a quadratic cost
diag_Q1 = @SVector [ # Player 1 state cost
    1., 1., 1., 1.,
    0., 0., 0., 0.]
diag_Q2 = @SVector [ # Player 2 state cost
    0., 0., 0., 0.,
    1., 1., 1., 1.]
Q = [0.1*Diagonal(diag_Q1), # Players state costs
     0.1*Diagonal(diag_Q2)]
Qf = [1.0*Diagonal(diag_Q1),
      1.0*Diagonal(diag_Q2)]

# Players controls costs
R = [0.1*Diagonal(@SVector ones(length(pu[1]))),
     0.1*Diagonal(@SVector ones(length(pu[2])))]

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
actors_types = [:car, :car]

# Create constraints
conSet = ConstraintSet(n,m,N)
con_inds = 2:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(conSet, actors_radii, px, p, con_inds)

# Define the problem
prob = GameProblem(model, obj, conSet, x0, xf, Z, N, tf)
# Specify the solver's options
opts = DirectGamesSolverOptions{T}()

solver = DirectGamesSolver(prob, opts)
reset!(solver, reset_type=:full)
solve!(solver)
```

## Examples
Notebooks with more detailed examples can be found [here](https://github.com/RoboticExplorationLab/ALGAMES.jl/tree/master/experiments/notebooks), including all the examples from our [RSS 2020 paper](https://github.com/RoboticExplorationLab/ALGAMES.jl/tree/master/experiments/rss_2020). Among these examples we solve a ramp merging problem and an intersection crossing problem:
![](plots/gif/intersection_3pl_2x.gif)
![](plots/gif/merging_3pl_2x.gif)

We also perform a Monte Carlo analysis of ALGAMES to evaluate the robustness and speed of the solver to noise in the initial condition.
![](plots/gif/monte_carlo_merging.gif)


<!-- ## Documentation
Detailed documentation for getting started with the package can be found [here](https://roboticexplorationlab.github.io/TrajectoryOptimization.jl/dev/). -->
