using Test

solver_2_players = GameProblems.ilqgames_ramp_merging_2_players_solver
solver_3_players = GameProblems.ilqgames_ramp_merging_3_players_solver
solver_4_players = GameProblems.ilqgames_ramp_merging_4_players_solver

# Solve with ALGAMES
solve!(solver_2_players)
solve!(solver_3_players)
solve!(solver_4_players)
@test converged(solver_2_players)
@test converged(solver_3_players)
@test converged(solver_4_players)
