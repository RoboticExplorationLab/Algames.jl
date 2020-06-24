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

Random.seed!(100)
xf_noise = 0.0*SVector{n}([# p1   # p2   # p3   # p4
					  0.000,  0.000,  0.000,  0.000, # x
					  0.010,  0.010,  0.010,  0.010, # y
					  0.010,  0.010,  0.010,  0.010, # θ
					  0.005,  0.005,  0.005,  0.005, # v
		            ])
# xfs = [SVector{n}(xf .+ xf_noise.*(rand(n) .- 0.5)) for i=1:p]
xfs = [SVector{n}(xf) for i=1:p]
# Players objectives
obj = [[LQRObjective(Q[j],R[j],Qf[j],xfs[i],N,checks=false) for j=1:p] for i=1:p]

# Create constraints
algames_conSet = [ConstraintSet(n,m,N) for i=1:p]
con_inds = 1:N # Indices where the constraints will be applied

for i = 1:p
	# # Add collision avoidance constraints only for the ego vehicle.
	# for j = 1:p
	# 	if j != i
	# 		add_collision_avoidance(algames_conSet[i], actors_radii, i, j, px,
	# 			p, con_inds; constraint_type=:constraint)
	# 	end
	# end
	# Add collision avoidance constraints
	add_collision_avoidance(algames_conSet[i], actors_radii, px,
		p, con_inds; constraint_type=:constraint)
	# Add scenario specific constraints
	add_scenario_constraints(algames_conSet[i], ramp_merging_4_players_unicycle_penalty_scenario,
		px, con_inds; constraint_type=:constraint)

	# Add controls constraints
	# u_lim = 0.12*ones(m)
	u_lim = U_lim*ones(m)
	control_bound = BoundConstraint(n,m,u_min=u_min,u_max=u_max)
	add_constraint!(algames_conSet[i], control_bound, 1:N-1)
end

approx_probs = [GameProblem(model, obj[i], xf, tf,
	constraints=algames_conSet[i], x0=x0, N=N) for i=1:p]

algames_opts = DirectGamesSolverOptions{T}(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
	optimality_constraint_tolerance=1e-2,
	μ_penalty=1.0,
    log_level=TO.Logging.Debug)
approx_solvers = [DirectGamesSolver(approx_probs[i], algames_opts) for i=1:p]

for i = 1:p
	# add penalty constraints
	add_collision_avoidance(approx_solvers[i].penalty_constraints,
	    inflated_actors_radii, px, p, con_inds; constraint_type=:constraint)
	reset!(approx_solvers[i], reset_type=:full)
end
@time timing_solve.(approx_solvers)
visualize_trajectory_car.(approx_solvers)
visualize_control.(approx_solvers)

# visualize_trajectory_car.(approx_solvers)
# #
# using MeshCat
# vis = MeshCat.Visualizer()
# anim = MeshCat.Animation()
# open(vis)
# sleep(1.0)
# # # Execute this line after the MeshCat tab is open
# vis, anim = animation(approx_solvers[1],
# 	ramp_merging_4_players_unicycle_penalty_scenario;
# 	vis=vis, anim=anim,
# 	open_vis=false,
# 	display_actors=true,
# 	display_trajectory=true)
