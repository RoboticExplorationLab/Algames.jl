using Test
using MeshCat

prob = GameProblems.algames_straight_2_players_prob
scenario = GameProblems.straight_2_players_scenario

opts = DirectGamesSolverOptions(
    iterations=10,
    inner_iterations=20,
    iterations_linesearch=10,
    min_steps_per_iteration=1,
    log_level=TO.Logging.Warn)
solver = DirectGamesSolver(prob, opts)
solve!(solver)

vis = MeshCat.Visualizer()
anim = MeshCat.Animation()
# open(vis)
# Execute this line after the MeshCat tab is open
vis, anim = animation(solver, scenario;
	vis=vis, anim=anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true)

@test converged(solver)
