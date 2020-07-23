export
    Scenario,
    IntersectionScenario,
    StraightScenario,
    TIntersectionScenario,
    MergingScenario


abstract type Scenario{T}
end
abstract type AbstractMergingScenario{T} <: Scenario{T}
end

@with_kw mutable struct IntersectionScenario{T} <:Scenario{T}
    l1::T # width of the roadway
    l2::T # length of the roadway
    actors_radii::Vector{T} # radii of the actors
    actors_types::Vector{Symbol} # type of the actors
    bound_radius::T # radius of the boundary corners
end

@with_kw mutable struct StraightScenario{T} <:Scenario{T}
    road_length::T=3.0 # width of the roadway
    road_width::T=0.37 # length of the roadway
    actors_radii::Vector{T} # radii of the actors
    actors_types::Vector{Symbol} # type of the actors
    obs_radius::T=0.06 # radius of the obstacle
    obs::Vector{T}=[0.50, -0.17] # 2D posiion of the obstacle
end

@with_kw mutable struct TIntersectionScenario{T} <: Scenario{T}
    top_road_length::T=4.0 # length of the top roadway
    road_width::T=0.60 # width of the roadway
    bottom_road_length::T=1.0 # length of the bottom roadway
    cross_width::T=0.25 # length of the bottom roadway
    actors_radii::Vector{T} # radii of the actors (cars, pedestrians)
    actors_types::Vector{Symbol} # type of the actors
    bound_radius::T=0.05 # radius of the boundary conrners
end

@with_kw mutable struct MergingScenario{T} <: AbstractMergingScenario{T}
    road_length::T=6.0 # length of the roadway
    road_width::T=0.34 # width of the roadway
    ramp_length::T=3.2 # length of the ramp
    ramp_angle::T=pi/12 # angle of the ramp wrt to the min roadway
    actors_radii::Vector{T} # radii of the actors (cars)
    actors_types::Vector{Symbol} # type of the actors
end

@with_kw mutable struct MergingObstacleScenario{T} <: AbstractMergingScenario{T}
    road_length::T=6.0 # length of the roadway
    road_width::T=0.34 # width of the roadway
    ramp_length::T=3.2 # length of the ramp
    ramp_angle::T=pi/12 # angle of the ramp wrt to the min roadway
    actors_radii::Vector{T} # radii of the actors (cars)
    actors_types::Vector{Symbol} # type of the actors
    obstacles_centers::Vector{Vector{T}} # center of the obstacles
    obstacles_radii::Vector{T} # radii of the obstacles
end

function Base.copy(s::S) where {S<:Scenario}
    fnames = fieldnames(typeof(s))
    args = [getfield(s,fname) for fname in fnames]
    S(args...)
end
