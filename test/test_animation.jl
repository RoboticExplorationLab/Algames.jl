using Test
using MeshCat

straight_prob = GameProblems.algames_straight_2_players_prob
ramp_merging_prob = GameProblems.algames_ramp_merging_2_players_prob
t_intersection_prob = GameProblems.algames_t_intersection_2_players_prob

straight_scenario = GameProblems.straight_2_players_scenario
ramp_merging_scenario = GameProblems.ramp_merging_2_players_scenario
t_intersection_scenario = GameProblems.t_intersection_2_players_scenario

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)

straight_solver = DirectGamesSolver(straight_prob, opts)
ramp_merging_solver = DirectGamesSolver(ramp_merging_prob, opts)
t_intersection_solver = DirectGamesSolver(t_intersection_prob, opts)

solve!(straight_solver)
solve!(ramp_merging_solver)
solve!(t_intersection_solver)

vis = MeshCat.Visualizer()
anim = MeshCat.Animation()
# open(vis)
# Execute this line after the MeshCat tab is open
vis, anim = animation(straight_solver, straight_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

vis, anim = animation(ramp_merging_solver, ramp_merging_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

vis, anim = animation(t_intersection_solver, t_intersection_scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

@test converged(straight_solver)
@test converged(ramp_merging_solver)
@test converged(t_intersection_solver)
