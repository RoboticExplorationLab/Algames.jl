export
	build_scenario,
	add_scenario_constraints,
	add_t_intersection_constraints,
	t_intersection_landmarks,
	add_rounded_boundary_constraints

function build_scenario(vis::Visualizer, scenario::TIntersectionScenario{T}; scale::T=10.) where T
	pkg_path = joinpath(dirname(@__FILE__), "../../../")
    # Plot Road in Meshcat
    road_image = PngImage(joinpath(pkg_path, "resources/textures/road.png"))
    road_texture = Texture(image=road_image)
    road_material = MeshLambertMaterial(map=road_texture)
	thickness = 0.002*scale

	road_width = scenario.road_width*scale
	top_length = scenario.top_road_length*scale
	bot_length = scenario.bottom_road_length*scale
	cross_width = scenario.cross_width*scale

	top_road = HyperRectangle(Vec(-top_length/2, -road_width/2, -thickness),
		Vec(top_length, road_width, thickness))
	bot_road = HyperRectangle(Vec(-road_width/2, -road_width/2-bot_length, -thickness),
		Vec(road_width, bot_length, thickness))
	setobject!(vis["roadway/top_road"], top_road, road_material)
	setobject!(vis["roadway/bot_road"], bot_road, road_material)

    # Plot lines in Meshcat
    line_material = MeshPhongMaterial(color=RGBA(1, 1, 0, 1.0))
	line_width = 0.005*scale
	top_left_line = HyperRectangle(Vec(-top_length/2, -line_width/2, 0.),
		Vec((top_length-road_width)/2, line_width, thickness))
	top_right_line = HyperRectangle(Vec(road_width/2, -line_width/2, 0.),
		Vec((top_length-road_width)/2, line_width, thickness))
	bottom_line = HyperRectangle(Vec(-line_width/2, -road_width/2-bot_length, 0.),
		Vec(line_width, bot_length, thickness))
	setobject!(vis["roadway/top_left_line"], top_left_line, line_material)
	setobject!(vis["roadway/top_right_line"], top_right_line, line_material)
	setobject!(vis["roadway/bottom_line"], bottom_line, line_material)

	# PLot boundaries in MeshCat
	bound_width = 0.015*scale
	bound_height = 0.03*scale
	bound_image = PngImage(joinpath(pkg_path, "resources/textures/light_boundary.png"))
	bound_texture = Texture(image=bound_image)
	bound_material = MeshLambertMaterial(map=bound_texture)
	# Upper Boundary
	up_bound = HyperRectangle(Vec(-top_length/2, road_width/2, 0.0),
		Vec(top_length, bound_width, bound_height))
	left_bound = HyperRectangle(Vec(-top_length/2, -road_width/2-bound_width, 0.0),
		Vec(top_length/2-road_width/2, bound_width, bound_height))
	right_bound = HyperRectangle(Vec(road_width/2, -road_width/2-bound_width, 0.0),
		Vec(top_length/2-road_width/2, bound_width, bound_height))
	bot_left_bound = HyperRectangle(Vec(-road_width/2-bound_width, -road_width/2-bot_length, 0.0),
		Vec(bound_width, bot_length, bound_height))
	bot_right_bound = HyperRectangle(Vec(road_width/2, -road_width/2-bot_length, 0.0),
		Vec(bound_width, bot_length, bound_height))

	setobject!(vis["roadway/up_bound"], up_bound, bound_material)
	setobject!(vis["roadway/left_bound"], left_bound, bound_material)
	setobject!(vis["roadway/right_bound"], right_bound, bound_material)
	setobject!(vis["roadway/bot_left_bound"], bot_left_bound, bound_material)
	setobject!(vis["roadway/bot_right_bound"], bot_right_bound, bound_material)

	# Plot crosswalks in MeshCat
	cross_height = 0.003*scale
	cross_image = PngImage(joinpath(pkg_path, "resources/textures/crosswalk.png"))
	cross_texture = Texture(image=cross_image)
	cross_material = MeshLambertMaterial(map=cross_texture)
	left_cross = HyperRectangle(Vec(-cross_width-road_width/2, -road_width/2, 0.0),
		Vec(cross_width, road_width, cross_height))
	right_cross = HyperRectangle(Vec(road_width/2, -road_width/2, 0.0),
		Vec(cross_width, road_width, cross_height))
	setobject!(vis["roadway/left_cross"], left_cross, cross_material)
	setobject!(vis["roadway/right_cross"], right_cross, cross_material)

	return nothing
end


# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_scenario_constraints(conSet::ConstraintSet,
	scenario::TIntersectionScenario{T}, lanes, px, con_inds; constraint_type=:constraint) where T
    for i = 1:length(px)
        add_scenario_constraints(conSet::ConstraintSet, scenario::TIntersectionScenario,
            lanes[i], i, px, con_inds; constraint_type=constraint_type)
    end
    return nothing
end


# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_scenario_constraints(conSet::ConstraintSet, scenario::TIntersectionScenario,
	lane, player_id, px, con_inds; constraint_type=:constraint)

	top_length = scenario.top_road_length
	road_width = scenario.road_width
	bot_length = scenario.bottom_road_length
	bound_radius = scenario.bound_radius
	P = t_intersection_landmarks(scenario, scenario.actors_radii[player_id]) #### take the largest radius
	V = [[1., 0.], [0., 1.]]
	if lane == 1
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[9]),
			SVector{2}(P[12]),
			SVector{2}(V[2]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[13], P[14],  P[19], bound_radius,
            player_id, px, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
			P[20], P[15],  P[18], bound_radius,
			player_id, px, con_inds)
	elseif lane == 2
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[1]),
			SVector{2}(P[2]),
			SVector{2}(V[2]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[3],  P[4],  P[19], bound_radius,
            player_id, px, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
			P[20], P[5],  P[8], bound_radius,
			player_id, px, con_inds)
	elseif lane == 3
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[1]),
			SVector{2}(P[2]),
			SVector{2}(V[2]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[3],  P[6],  P[21], bound_radius,
            player_id, px, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
			P[26], P[11],  P[12], bound_radius,
			player_id, px, con_inds)
		add_rounded_boundary_constraints(conSet::ConstraintSet,
			P[18], P[17],  P[22], bound_radius,
			player_id, px, con_inds)
	elseif lane == 4
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[23]),
			SVector{2}(P[29]),
			SVector{2}(-V[1]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[24]),
			SVector{2}(P[30]),
			SVector{2}(V[1]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
	elseif lane == 5
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[27]),
			SVector{2}(P[31]),
			SVector{2}(-V[1]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
		con = BoundaryConstraint(conSet.n,
			SVector{2}(P[28]),
			SVector{2}(P[32]),
			SVector{2}(V[1]),
			px[player_id]...)
		add_constraint!(conSet, con, con_inds)
	end
    return nothing
end


function t_intersection_landmarks(scenario::TIntersectionScenario, actor_radius)
	r = actor_radius
	top_length = scenario.top_road_length
	road_width = scenario.road_width
	bot_length = scenario.bottom_road_length
	cross_width = scenario.cross_width

	l6 = road_width/2
	l5 = road_width/2 + cross_width - r
	l4 = road_width/2 + r
	l3 = road_width/2 + bot_length
	l2 = top_length/2
	l1 = road_width/2 - r
	l0 = r
    P = [[-l2, l1], #Line1
         [ l2, l1],

         [-l2, l0],
		 [-l1, l0],
		 [-l0, l0],
		 [ l0, l0],
		 [ l1, l0],
		 [ l2, l0],

         [-l2,-l0], #Line2
		 [-l1,-l0],
		 [ l1,-l0],
		 [ l2,-l0],

         [-l2,-l1],
		 [-l1,-l1],
		 [-l0,-l1],
		 [ l0,-l1],
         [ l1,-l1],
         [ l2,-l1],

         [-l1,-l3], #Line3
         [-l0,-l3],
         [ l0,-l3],
         [ l1,-l3],

         [-l5, l6], #Line4
		 [-l4, l6],
		 [-l1, l6],
		 [ l1, l6],
		 [ l4, l6],
		 [ l5, l6],

         [-l5,-l6], #Line5
         [-l4,-l6],
         [ l4,-l6],
         [ l5,-l6]]
    return P
end

function add_rounded_boundary_constraints(conSet::ConstraintSet, p1, p2, p3, r,
    player_id, px, con_inds)
    # p3 *          y
    #    |          |
    #    |___.      ---> x
    #    |⋱ |r
    # p2 *----------------* p1

    x = (p1 - p2)./norm(p1 - p2)
    y = (p3 - p2)./norm(p3 - p2)
    # Bound 1
    b1 = SVector{2}(p1)
    b2 = p2 + r*x
    b2 = SVector{2}(b2)
    v = SVector{2}(y)
    con = BoundaryConstraint(conSet.n, b1, b2, v, px[player_id]...)
    add_constraint!(conSet, con, con_inds)
    # Bound 2
    b1 = p2 + r*y
    b1 = SVector{2}(b1)
    b2 = SVector{2}(p3)
    v = SVector{2}(x)
    con = BoundaryConstraint(conSet.n, b1, b2, v, px[player_id]...)
    add_constraint!(conSet, con, con_inds)
    # Circle
    c = p2 + r*(x+y)
    x = SVector{1}([c[1]])
    y = SVector{1}([c[2]])
    radii = SVector{1}([r])
    con = CircleConstraint(conSet.n, x, y, radii, px[player_id]...)
    add_constraint!(conSet, con, con_inds)
    return nothing
end
