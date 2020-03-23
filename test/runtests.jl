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

# include("../src/game_problems.jl")

@testset "algames" begin
    include("test_algames.jl")
end

@testset "ilqgames" begin
    include("test_ilqgames.jl")
end

@testset "mpc" begin
    include("test_mpc.jl")
end

@testset "sampler" begin
    include("test_sampler.jl")
end

@testset "lq_game" begin
    include("test_lq_game.jl")
end
