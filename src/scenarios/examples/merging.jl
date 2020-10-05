export
	build_scenario,
	add_scenario_constraints,
	add_merging_constraints

function build_scenario(vis::Visualizer, scenario::AbstractMergingScenario; scale::T=10.) where T
	pkg_path = joinpath(dirname(@__FILE__), "../../../")
    # Plot Road in Meshcat
    road_image = PngImage(joinpath(pkg_path, "resources/textures/road.png"))
    road_texture = Texture(image=road_image)
    road_material = MeshLambertMaterial(map=road_texture)
	thickness = 0.002*scale

	road_width = scenario.road_width*scale
	road_length = scenario.road_length*scale
	ramp_length = scenario.ramp_length*scale
	ramp_angle = scenario.ramp_angle
	ext_ramp_length = ramp_length / cos(ramp_angle)
	ramp_delta = ext_ramp_length * sin(ramp_angle)

	road = HyperRectangle(Vec(-road_length/2, -road_width/2, -thickness),
		Vec(road_length, road_width, thickness))
	ramp = HyperRectangle(Vec(0., 0., -thickness*1.5),
		Vec(ext_ramp_length, road_width, thickness))
	setobject!(vis["roadway/road"], road, road_material)
	setobject!(vis["roadway/ramp"], ramp, road_material)

	ramp_rotation = LinearMap(AngleAxis(ramp_angle, 0, 0, 1))
	ramp_translation = Translation(-road_length/2, -road_width/2-ramp_delta, 0.)
	ramp_transformation = compose(ramp_translation, ramp_rotation)
	settransform!(vis["roadway/ramp"], ramp_transformation)

    # Plot lines in Meshcat
    line_material = MeshPhongMaterial(color=RGBA(1, 1, 0, 1.0))
	line_width = 0.005*scale
	top_line = HyperRectangle(Vec(-road_length/2, -line_width/2, 0.),
		Vec(road_length, line_width, thickness))
	bot_line = HyperRectangle(Vec(-road_length/2, -line_width/2-road_width/2, 0.),
		Vec(ramp_length, line_width, thickness))
		setobject!(vis["roadway/top_line"], top_line, line_material)
		setobject!(vis["roadway/bot_line"], bot_line, line_material)

	# Plot boundaries in MeshCat
	bound_width = 0.015*scale
	bound_height = 0.03*scale
	bound_image = PngImage(joinpath(pkg_path, "resources/textures/light_boundary.png"))
	bound_texture = Texture(image=bound_image)
	bound_material = MeshLambertMaterial(map=bound_texture)
	obs_image = PngImage(joinpath(pkg_path, "resources/textures/obstacle.png"))
	obs_texture = Texture(image=obs_image)
	obs_material = MeshLambertMaterial(map=obs_texture)
	# Upper Boundary
	up_bound = HyperRectangle(Vec(-road_length/2, road_width/2, 0.),
		Vec(road_length, bound_width, bound_height))
	bot_bound = HyperRectangle(Vec(-road_length/2+ramp_length, -road_width/2-bound_width, 0.0),
		Vec(road_length-ramp_length, bound_width, bound_height))
	ramp_bound = HyperRectangle(Vec(0., -bound_width, 0.),
		Vec(ext_ramp_length, bound_width, bound_height))
	setobject!(vis["roadway/up_bound"], up_bound, bound_material)
	setobject!(vis["roadway/bot_bound"], bot_bound, bound_material)
	setobject!(vis["roadway/ramp_bound"], ramp_bound, bound_material)
	settransform!(vis["roadway/ramp_bound"], ramp_transformation)

	# Plot obstacles in MeshCat.
	obstacle_height = 1.25*bound_height
	if typeof(scenario) <: MergingObstacleScenario
		obs_centers = scenario.obstacles_centers
		obs_radii = scenario.obstacles_radii
		s = length(obs_radii)
		for k = 1:s
			# Add obstacle cylinders
			obstacle_cylinder = Cylinder(
				Point(scale*obs_centers[k][1],scale*obs_centers[k][2],0.0),
				Point(scale*obs_centers[k][1],scale*obs_centers[k][2],obstacle_height),
				scale*obs_radii[k])
			setobject!(vis["roadway/obs_cyl_$k"], obstacle_cylinder, obs_material)
		end
	end
	return nothing
end

# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_scenario_constraints(conSet::ConstraintSet,
	scenario::AbstractMergingScenario, px, con_inds; constraint_type=:constraint)
	p = length(px)
    for i = 1:p
        add_scenario_constraints(conSet::ConstraintSet, scenario::AbstractMergingScenario,
            i, px, con_inds; constraint_type=constraint_type)
    end
	# Obstacle constraints
	if typeof(scenario) <: MergingObstacleScenario
		add_obstacle_avoidance!(
			conSet,
			scenario.obstacles_centers,
			scenario.obstacles_radii,
			scenario.actors_radii,
			px,
			con_inds)
	end
    return nothing
end


# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_scenario_constraints(conSet::ConstraintSet, scenario::AbstractMergingScenario,
    player_id, px, con_inds; constraint_type=:constraint)
	road_length = scenario.road_length
	road_width = scenario.road_width
	ramp_length = scenario.ramp_length
	ramp_angle = scenario.ramp_angle
	ext_ramp_length = ramp_length / cos(ramp_angle)
	ramp_delta = ext_ramp_length * sin(ramp_angle)

    l0 = road_width/2 - scenario.actors_radii[player_id] ### only take into account the last radius
    l2 = road_length/2
    b1 = @SVector [-l2,  l0]
    b2 = @SVector [ l2,  l0]
    b3 = @SVector [ l2, -l0]
	b4 = @SVector [-l2+ramp_length, -l0]
	b5 = @SVector [-l2, -l0-ramp_delta]
    v1 = @SVector [0.,  1.]
	v2 = @SVector [0., -1.]
	v3 = @SVector [sin(ramp_angle), -cos(ramp_angle)]

    con = BoundaryConstraint(conSet.n, b1, b2, v1, px[player_id]...)
    add_constraint!(conSet, con, con_inds)
    con = BoundaryConstraint(conSet.n, b3, b4, v2, px[player_id]...)
    add_constraint!(conSet, con, con_inds)
	con = BoundaryConstraint(conSet.n, b4, b5, v3, px[player_id]...)
    add_constraint!(conSet, con, con_inds)

    return nothing
end
