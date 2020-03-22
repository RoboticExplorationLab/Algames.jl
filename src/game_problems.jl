"""
    Game Problems
        Collection of game problems
"""
module GameProblems

using ALGAMES
using LinearAlgebra
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
const AG = ALGAMES

# include("../game_problems/ramp_merging.jl")
include("../game_problems/ramp_merging_2_players.jl")
include("../game_problems/ramp_merging_3_players.jl")
include("../game_problems/ramp_merging_4_players.jl")
include("../game_problems/t_intersection_2_players.jl")
include("../game_problems/t_intersection_3_players.jl")
include("../game_problems/t_intersection_4_players.jl")

export
    algames_ramp_merging_2_players_prob,
    algames_ramp_merging_3_players_prob,
    algames_ramp_merging_4_players_prob,
    ilqgames_ramp_merging_2_players_prob,
    ilqgames_ramp_merging_3_players_prob,
    ilqgames_ramp_merging_4_players_prob

export
    algames_ramp_merging_2_players_solver,
    algames_ramp_merging_3_players_solver,
    algames_ramp_merging_4_players_solver,
    ilqgames_ramp_merging_2_players_solver,
    ilqgames_ramp_merging_3_players_solver,
    ilqgames_ramp_merging_4_players_solver

export
    algames_t_intersection_2_players_prob,
    algames_t_intersection_3_players_prob,
    algames_t_intersection_4_players_prob,
    ilqgames_t_intersection_2_players_prob,
    ilqgames_t_intersection_3_players_prob,
    ilqgames_t_intersection_4_players_prob

export
    algames_t_intersection_2_players_solver,
    algames_t_intersection_3_players_solver,
    algames_t_intersection_4_players_solver,
    ilqgames_t_intersection_2_players_solver,
    ilqgames_t_intersection_3_players_solver,
    ilqgames_t_intersection_4_players_solver

end
