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

include("../game_problems/ramp_merging.jl")

export
    algames_ramp_merging,
    ilqqgames_ramp_merging

end
