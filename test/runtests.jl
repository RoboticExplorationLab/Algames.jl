using ALGAMES
using TrajectoryOptimization
using Test
using BenchmarkTools
using LinearAlgebra
using Random
using StaticArrays
using SparseArrays
using Logging
const TO = TrajectoryOptimization
const AG = ALGAMES

include("../src/game_problems.jl")

@testset "algames" begin
    include("test_algames.jl")
end

@testset "ilqgames" begin
    include("test_ilqgames.jl")
end
