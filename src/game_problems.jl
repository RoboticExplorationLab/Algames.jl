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

include("../game_problems/linear_quadratic.jl")

include("../game_problems/mpc_ramp_merging_2_players.jl")
include("../game_problems/mpc_ramp_merging_3_players.jl")

include("../game_problems/lane_changing_3_players_penalty.jl")
include("../game_problems/overtaking_3_players_penalty.jl")

include("../game_problems/ramp_merging_2_players_penalty.jl")
include("../game_problems/ramp_merging_3_players_penalty.jl")
include("../game_problems/ramp_merging_4_players_penalty.jl")

include("../game_problems/ramp_merging_2_players_unicycle_penalty.jl")
include("../game_problems/ramp_merging_3_players_unicycle_penalty.jl")
include("../game_problems/ramp_merging_4_players_unicycle_penalty.jl")

include("../game_problems/t_intersection_3_players_penalty.jl")

include("../game_problems/ramp_merging_2_players.jl")
include("../game_problems/ramp_merging_3_players.jl")
include("../game_problems/ramp_merging_4_players.jl")

include("../game_problems/straight_2_players.jl")

include("../game_problems/t_intersection_2_players.jl")
include("../game_problems/t_intersection_3_players.jl")
include("../game_problems/t_intersection_4_players.jl")

export
    algames_linear_quadratic_prob,
    ilqgames_linear_quadratic_prob,

export
    ramp_merging_2_players_mpc_scenario,
    ramp_merging_3_players_mpc_scenario,
    algames_ramp_merging_2_players_mpc_prob,
    algames_ramp_merging_3_players_mpc_prob,
    algames_ramp_merging_2_players_mpc_opts,
    algames_ramp_merging_3_players_mpc_opts,
    ramp_merging_2_players_mpc_opts,
    ramp_merging_3_players_mpc_opts

export
    lane_changing_3_players_penalty_scenario,
    algames_lane_changing_3_players_penalty_prob,
    algames_lane_changing_3_players_penalty_opts,
    algames_lane_changing_3_players_penalty_contraints

export
    overtaking_3_players_penalty_scenario,
    algames_overtaking_3_players_penalty_prob,
    algames_overtaking_3_players_penalty_opts,
    algames_overtaking_3_players_penalty_contraints

export
    ramp_merging_2_players_penalty_scenario,
    ramp_merging_3_players_penalty_scenario,
    ramp_merging_4_players_penalty_scenario,
    algames_ramp_merging_2_players_penalty_prob,
    algames_ramp_merging_3_players_penalty_prob,
    algames_ramp_merging_4_players_penalty_prob,
    algames_ramp_merging_2_players_penalty_opts,
    algames_ramp_merging_3_players_penalty_opts,
    algames_ramp_merging_4_players_penalty_opts,
    algames_ramp_merging_2_players_penalty_contraints,
    algames_ramp_merging_3_players_penalty_contraints,
    algames_ramp_merging_4_players_penalty_contraints


export
    ramp_merging_2_players_unicycle_penalty_scenario,
    ramp_merging_3_players_unicycle_penalty_scenario,
    ramp_merging_4_players_unicycle_penalty_scenario,
    algames_ramp_merging_2_players_unicycle_penalty_prob,
    algames_ramp_merging_3_players_unicycle_penalty_prob,
    algames_ramp_merging_4_players_unicycle_penalty_prob,
    algames_ramp_merging_2_players_unicycle_penalty_opts,
    algames_ramp_merging_3_players_unicycle_penalty_opts,
    algames_ramp_merging_4_players_unicycle_penalty_opts,
    algames_ramp_merging_2_players_unicycle_penalty_contraints,
    algames_ramp_merging_3_players_unicycle_penalty_contraints,
    algames_ramp_merging_4_players_unicycle_penalty_contraints

export
    t_intersection_3_players_penalty_scenario,
    algames_t_intersection_3_players_penalty_prob,
    algames_t_intersection_3_players_penalty_opts,
    algames_t_intersection_3_players_penalty_contraints

export
    straight_2_players_scenario,
    algames_straight_2_players_prob,
    ilqgames_straight_2_players_prob

export
    ramp_merging_2_players_scenario,
    ramp_merging_3_players_scenario,
    ramp_merging_4_players_scenario,
    algames_ramp_merging_2_players_prob,
    algames_ramp_merging_3_players_prob,
    algames_ramp_merging_4_players_prob,
    ilqgames_ramp_merging_2_players_prob,
    ilqgames_ramp_merging_3_players_prob,
    ilqgames_ramp_merging_4_players_prob,
    algames_ramp_merging_2_players_opts,
    algames_ramp_merging_3_players_opts,
    algames_ramp_merging_4_players_opts,
    ilqgames_ramp_merging_2_players_opts,
    ilqgames_ramp_merging_3_players_opts,
    ilqgames_ramp_merging_4_players_opts

export
    t_intersection_2_players_scenario,
    t_intersection_3_players_scenario,
    t_intersection_4_players_scenario,
    algames_t_intersection_2_players_prob,
    algames_t_intersection_3_players_prob,
    algames_t_intersection_4_players_prob,
    ilqgames_t_intersection_2_players_prob,
    ilqgames_t_intersection_3_players_prob,
    ilqgames_t_intersection_4_players_prob,
    algames_t_intersection_2_players_opts,
    algames_t_intersection_3_players_opts,
    algames_t_intersection_4_players_opts,
    ilqgames_t_intersection_2_players_opts,
    ilqgames_t_intersection_3_players_opts,
    ilqgames_t_intersection_4_players_opts

end
