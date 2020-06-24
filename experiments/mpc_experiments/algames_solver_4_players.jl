using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

include("solver_4_players.jl")

Q  = [0.1*Diagonal(diag_Q[i]) for i=1:p] # Players stage state costs
Qf = [1.0*Diagonal(diag_Q[i]) for i=1:p] # Players final state costs
# Players controls costs
R = [0.1*Diagonal(diag_R[i]) for i=1:p]

# Players objectives
obj = [LQRObjective(Q[i],R[i],Qf[i],xf,N,checks=false) for i=1:p]

# Create constraints
algames_conSet = ConstraintSet(n,m,N)
con_inds = 1:N # Indices where the constraints will be applied

# Add collision avoidance constraints
add_collision_avoidance(algames_conSet, actors_radii, px,
	p, con_inds; constraint_type=:constraint)
# Add collision avoidance constraints########################
add_collision_avoidance(algames_conSet, 1.2*actors_radii, px,########################
	p, 7:N; constraint_type=:constraint)########################
# Add collision avoidance constraints########################
add_collision_avoidance(algames_conSet, 1.3*actors_radii, px,########################
	p, 14:N; constraint_type=:constraint)########################
# Add scenario specific constraints
add_scenario_constraints(algames_conSet, ramp_merging_4_players_unicycle_penalty_scenario,
	px, con_inds; constraint_type=:constraint)

# Add controls constraints
control_bound = BoundConstraint(n,m,u_min=u_min,u_max=u_max)
add_constraint!(algames_conSet, control_bound, 1:N-1)

algames_prob = GameProblem(model, obj, xf, tf,
	constraints=algames_conSet, x0=x0, N=N)

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
	optimality_constraint_tolerance=1e-2,
	μ_penalty=μ_penalty,
    log_level=TO.Logging.Debug)
algames_solver = DirectGamesSolver(algames_prob, algames_opts)

# add penalty constraints
add_collision_avoidance(algames_solver.penalty_constraints,
    inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)

reset!(algames_solver, reset_type=:full)
# algames_ramp_merging_2_players_unicycle_penalty_contraints = copy(algames_ramp_merging_2_players_unicycle_penalty_solver.penalty_constraints)

@time timing_solve(algames_solver)
visualize_trajectory_car(algames_solver)
# visualize_control(algames_solver)

# using MeshCat
# vis = MeshCat.Visualizer()
# anim = MeshCat.Animation()
# open(vis)
# sleep(1.0)
# # Execute this line after the MeshCat tab is open
# vis, anim = animation(algames_solver,
# 	ramp_merging_4_players_unicycle_penalty_scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=true)
