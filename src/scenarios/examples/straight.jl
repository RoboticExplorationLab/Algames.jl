export
	build_scenario,
	add_scenario_constraints,
	add_straight_constraints


function build_scenario(vis::Visualizer, scenario::StraightScenario{T}; scale::T=10.) where T
	pkg_path = joinpath(dirname(@__FILE__), "../../../")
    # Plot Road in Meshcat
    road_image = PngImage(joinpath(pkg_path, "resources/textures/road.png"))
    road_texture = Texture(image=road_image)
    road_material = MeshLambertMaterial(map=road_texture)
	thickness = 0.002*scale

	road_width = scenario.road_width*scale
	road_length = scenario.road_length*scale
	actors_radii = scenario.actors_radii*scale
	obs_radius = scenario.obs_radius*scale
	obs = scenario.obs*scale

	actor_radius = actors_radii[1]
	# Check all radii are the same
	check = all(actors_radii .== actor_radius) ### only check the last radius
	if !check
		@warn "All actors radii are not equal this
		will create inconsistent vizualizations."
	end

	road = HyperRectangle(Vec(-road_length/2, -road_width/2, -thickness),
		Vec(road_length, road_width, thickness))
	setobject!(vis["roadway/road"], road, road_material)

    # Plot lines in Meshcat
    line_material = MeshPhongMaterial(color=RGBA(1, 1, 0, 1.0))
	line_width = 0.005*scale
	line = HyperRectangle(Vec(-road_length/2, -line_width/2, 0.),
		Vec(road_length, line_width, thickness))
	setobject!(vis["roadway/line"], line, line_material)

	# Plot boundaries in MeshCat
	bound_width = 0.015*scale
	bound_height = 0.03*scale
	bound_image = PngImage(joinpath(pkg_path, "resources/textures/dark_boundary.png"))
	bound_texture = Texture(image=bound_image)
	bound_material = MeshLambertMaterial(map=bound_texture)

	up_bound = HyperRectangle(Vec(-road_length/2, road_width/2, 0.),
		Vec(road_length, bound_width, bound_height))
	bot_bound = HyperRectangle(Vec(-road_length/2, -road_width/2-bound_width, 0.0),
		Vec(road_length, bound_width, bound_height))
	setobject!(vis["roadway/up_bound"], up_bound, bound_material)
	setobject!(vis["roadway/bot_bound"], bot_bound, bound_material)

	# Plot cylindrical obstacle
	obs_height = 0.03*scale
	obstacle = Cylinder(Point(0., 0., 0.),
		Point(0., 0., obs_height), obs_radius)
	setobject!(vis["roadway/obs"], obstacle, bound_material)
	obs_translation = Translation(obs..., 0.)
	settransform!(vis["roadway/obs"], obs_translation)

	return nothing
end


# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_scenario_constraints(conSet::ConstraintSet, scenario::StraightScenario, px,
    con_inds; constraint_type=:constraint)
    for i = 1:length(px)
        add_scenario_constraints(conSet::ConstraintSet,
            scenario.road_length,
            scenario.road_width,
            scenario.actors_radii[i],
			scenario.obs,
			scenario.obs_radius,
            i, px, con_inds)
    end
    return nothing
end


# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_scenario_constraints(conSet::ConstraintSet, road_length, road_width, car_radius,
	obs, obs_radius, player_id, px, con_inds; constraint_type=:constraint)
    l0 = road_width/2 - car_radius
    l2 = road_length/2
    b1 = @SVector [-l2,  l0]
    b2 = @SVector [ l2,  l0]
    b3 = @SVector [ l2, -l0]
    b4 = @SVector [-l2, -l0]
    v1 = @SVector [0.,  1.]
    v2 = @SVector [0., -1.]
    con = BoundaryConstraint(conSet.n, b1, b2, v1, px[player_id]...)
    add_constraint!(conSet, con, con_inds)
    con = BoundaryConstraint(conSet.n, b3, b4, v2, px[player_id]...)
    add_constraint!(conSet, con, con_inds)

	if constraint_type == :constraint
		con = CircleConstraint(conSet.n,
			SVector{1}([obs[1]]),
			SVector{1}([obs[2]]),
			SVector{1}([obs_radius+car_radius]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
	elseif constraint_type == :penalty
		con = ExpCircleConstraint(conSet.n,
			SVector{1}([obs[1]]),
			SVector{1}([obs[2]]),
			SVector{1}([obs_radius+car_radius]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
	else
		@warn "Wrong type of constraint."
	end
    return nothing
end
