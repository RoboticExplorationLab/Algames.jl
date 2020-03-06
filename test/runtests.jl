using ALGAMES
using TrajectoryOptimization
using Test
using BenchmarkTools
using LinearAlgebra
using Random
using StaticArrays
using SparseArrays
# using Ipopt
using Logging
# using ForwardDiff
const TO = TrajectoryOptimization
const AG = ALGAMES

@testset "algames" begin
    include("test_algames.jl")
end
