export
    Scenario,
    IntersectionScenario,
    StraightScenario,
    TIntersectionScenario,
    MergingScenario


abstract type Scenario{T}
end

struct IntersectionScenario{T} <:Scenario{T}
    l1::T # width of the roadway
    l2::T # length of the roadway
    actors_radii::Vector{T} # radii of the actors
    actors_types::Vector{Symbol} # type of the actors
    bound_radius::T # radius of the boundary corners
end

struct StraightScenario{T} <:Scenario{T}
    road_length::T # width of the roadway
    road_width::T # length of the roadway
    actors_radii::Vector{T} # radii of the actors
    actors_types::Vector{Symbol} # type of the actors
    obs_radius::T # radius of the obstacle
    obs::Vector{T} # 2D posiion of the obstacle
end

struct TIntersectionScenario{T} <: Scenario{T}
    top_road_length::T # length of the top roadway
    road_width::T # width of the roadway
    bottom_road_length::T # length of the bottom roadway
    cross_width::T # length of the bottom roadway
    actors_radii::Vector{T} # radii of the actors (cars, pedestrians)
    actors_types::Vector{Symbol} # type of the actors
    bound_radius::T # radius of the boundary conrners
end

struct MergingScenario{T} <: Scenario{T}
    road_length::T # length of the roadway
    road_width::T # width of the roadway
    ramp_length::T # length of the ramp
    ramp_angle::T # angle of the ramp wrt to the min roadway
    actors_radii::Vector{T} # radii of the actors (cars)
    actors_types::Vector{Symbol} # type of the actors
end
