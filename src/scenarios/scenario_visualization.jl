export
	solver_scope,
	animation,
	build_actors,
	build_trajectory,
	scene_animation,
	actor_transformations

function solver_scope(solver::DirectGamesSolver)
	return :DirectGamesSolver
end

function solver_scope(solver::iLQGamesSolver)
	return :iLQGamesSolver
end

function solver_scope(solver::PenaltyiLQGamesSolver)
	return :PenaltyiLQGamesSolver
end

function solver_scope(solver::MPCGamesSolver)
	return :MPCGamesSolver
end

function animation(solver::TO.AbstractSolver, scenario::Scenario{T};
	vis::Visualizer=Visualizer(), anim::MeshCat.Animation=MeshCat.Animation(),
	open_vis::Bool=true,
	display_actors::Bool=true,
	display_trajectory::Bool=false,
	no_background::Bool=false) where T

	# vis.core.scope
	delete!(vis["/Grid"])
	delete!(vis["/Axes"])
	no_background ? delete!(vis["Background"]) : nothing
	open_vis ? open(vis) : nothing
	display_actors ?
		build_actors(vis, scenario, solver_scope(solver)) :
		delete!(vis.core.tree.children["meshcat"].children[string(solver_scope(solver))], "actors")
	display_trajectory ?
		build_trajectory(vis, solver) :
		delete!(vis.core.tree.children["meshcat"].children[string(solver_scope(solver))], "trajectory")
	build_scenario(vis, scenario)

	# Animate the scene
	default_framerate = 3
	anim = MeshCat.Animation(anim.clips, default_framerate)

	scene_animation(solver, vis, anim)
    MeshCat.setanimation!(vis, anim)
	return vis, anim
end

function build_actors(vis::Visualizer, scenario::Scenario{T}, scope::Symbol) where T
	actors_radii = scenario.actors_radii
	actors_types = scenario.actors_types
	car_offset = compose(compose(compose(
			Translation(0.0, 0.0, -0.0),
			LinearMap(AngleAxis(pi/200, 0, 1, 0))),
			LinearMap(AngleAxis(pi/2, 0, 0, 1))),
			LinearMap(AngleAxis(pi/2, 1, 0, 0)))
	ped_offset = LinearMap(AngleAxis(pi/2, 0, 0, 1))

	cylinder_height_car = 0.01
	cylinder_height_ped = 0.035
	color_values = [100/255 149/255 237/255; # cornflowerblue
					1       0       0;       # red
					1       1       0;       # yellow
					34/255  139/255 34/255;  # forestgreen
					1       165/255 0;       # orange
					0.5     0.5     0.5;     # gray
					1       165/255 0;       # orange
					0       1       0]       # Lime
	colors = ["cornflowerblue", "red", "yellow", "forestgreen",
		"orange", "gray", "orange", "lime"]

	path = String(scope)*"/actors/actor_bundle_"
	for i = 1:length(actors_radii)
		if actors_types[i] == :car
			actor = load("resources/object/car_geometry.obj", GLUVMesh)
			actor.vertices .*= 0.025
			actor_image = PngImage("resources/textures/"*colors[i]*".png")
			actor_texture = Texture(image=actor_image)
			actor_material = MeshLambertMaterial(map=actor_texture)
			setobject!(vis[path*"$i/actor"], actor, actor_material)
			settransform!(vis[path*"$i/actor"], car_offset)

			# Add collision avoidance cylinders
			collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,cylinder_height_car), actors_radii[i]) ### need to plaot cylinders per timestep
			cylinder_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., 9e-1))
			setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)

		elseif actors_types[i] == :pedestrian
			actor = load("resources/object/pedestrian_geometry.obj", GLUVMesh)
			actor.vertices .*= 0.001/15
			actor_image = PngImage("resources/textures/black_boundary.png")
			actor_texture = Texture(image=actor_image)
			actor_material = MeshLambertMaterial(map=actor_texture)
			setobject!(vis[path*"$i/actor"], actor, actor_material)
			settransform!(vis[path*"$i/actor"], ped_offset)

			# Add collision avoidance cylinders
			collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,cylinder_height_ped), actors_radii[i]) ### need to plaot cylinders per timestep
			cylinder_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., 6e-1))
			setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)
		end
	end
	return nothing
end

function build_trajectory(vis::Visualizer, solver::TO.AbstractSolver{T}) where {T}
	color_values = [100/255 149/255 237/255; # cornflowerblue
					1       0       0;       # red
					1       1       0;       # yellow
					34/255  139/255 34/255;  # forestgreen
					1       165/255 0;       # orange
					0.5     0.5     0.5;     # gray
					1       165/255 0;       # orange
					0       1       0]       # Lime
	colors = ["cornflowerblue", "red", "yellow", "forestgreen",
		"orange", "gray", "orange", "lime"]
	cyl_height = 0.005
	cyl_radius = 0.025

	n,m,N = size(solver)
	model = TO.get_model(solver)
	n,m,pu,p = size(model)
	px = model.px
	path = String(solver_scope(solver))*"/trajectory"
	X = TO.states(solver)
	for k = 1:length(X)
		for i = 1:p
			cyl_pos = Translation(X[k][px[i]]..., 0)
			aug_cyl_height = cyl_height + i*0.00 + k*0.07/length(X)
			cyl = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,aug_cyl_height), cyl_radius)
			cyl_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., 9e-1))
			setobject!(vis[path*"/cyl_$i"*"_$k"], cyl, cyl_material)
			settransform!(vis[path*"/cyl_$i"*"_$k"], cyl_pos)
		end
	end
	return nothing
end

function scene_animation(solver::TO.AbstractSolver, vis::Visualizer, anim::MeshCat.Animation)
	n,m,N = size(solver)
	model = TO.get_model(solver)
	n,m,pu,p = size(model)
	X = TO.states(solver)
	U = TO.controls(solver)

	len = length(solver.Z)

	# Animate the scene
	birdseye_trans = Translation(0.0, 0.0, 1.1)
	# birdseye_trans = Translation(0.0, 0.0, 2.2)
    birdseye_rot = compose(
		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))

    # Compute actor transformations
    actor_translations, actor_rotations = actor_transformations(solver, len)
	# We get rid of the last frame if the final state constraints are relaxed
	path = String(solver_scope(solver))*"/actors/actor_bundle"
	for k=1:len-1
        # Set the poses of the two actors.
        MeshCat.atframe(anim, k) do
            for i=1:p
                settransform!(vis[path*"_$i"],
					compose(actor_translations[k][i],actor_rotations[k][i]))
            end
			setprop!(vis["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
			setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 1.0)
            camera_transformation = compose(
                birdseye_trans,
                birdseye_rot)
            settransform!(vis["/Cameras/default"], camera_transformation)
        end
	end
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 5.0)
	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
	return nothing
end

function actor_transformations(solver, iter)
	n,m,N = size(solver)
	px = TO.get_model(solver).px
	p = length(px)
	X = TO.states(solver)
    actor_translations = []
    actor_rotations = []
    for k=1:iter
        actor_translation = [Translation(X[k][px[i]]..., 0) for i=1:p]
        if k==iter
            # Deals with the last step, we assume constant heading
            actor_rotation = actor_rotations[end]
        else
            angles = [atan(X[k+1][px[i][2]]-X[k][px[i][2]],
					  X[k+1][px[i][1]]-X[k][px[i][1]]) for i=1:p]
            actor_rotation = [LinearMap(AngleAxis(angles[i], 0, 0, 1.)) for i=1:p]
        end
        push!(actor_translations, actor_translation)
        push!(actor_rotations, actor_rotation)
    end
    return actor_translations, actor_rotations
end




################################
################################
################################


function animation(solver_leader::TO.ALTROSolver{T},
	solver_follower::TO.ALTROSolver{T},
	scenario::Scenario{T}, plx;
	no_background=false) where T

	n,m,N = size(solver_leader)
    # Open visualizer
    vis = Visualizer()
	vis.core.scope
    open(vis)
	build_actors(vis, scenario)
	build_scenario(vis, scenario)
	# Animate the scene
	vis, anim = scene_animation(solver_leader, solver_follower,
		vis, plx)
    MeshCat.setanimation!(vis, anim)
end

function animation(solver::TO.ALTROSolver{T}, scenario::Scenario{T}, plx;
	no_background=false) where T

	n,m,N = size(solver)
    # Open visualizer
    vis = Visualizer()
	vis.core.scope
    open(vis)
	build_actors(vis, scenario)
	build_scenario(vis, scenario)
	# Animate the scene
	vis, anim = scene_animation(solver, vis, plx; no_background=no_background)
    MeshCat.setanimation!(vis, anim)
end

function scene_animation(solver::TO.ALTROSolver{T}, vis, plx; no_background=false) where T
	n,m,N = size(solver)
	X = TO.states(solver)
	U = TO.controls(solver)

	# Animate the scene
	birdseye_trans = Translation(0.0, 0.0, -2.2)
    birdseye_rot = compose(
		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))

	# Plot Trajectory
    anim = MeshCat.Animation()
	default_framerate = 3
	anim = MeshCat.Animation(anim.clips, default_framerate)
	delete!(vis["/Grid"])
	delete!(vis["/Axes"])
	if no_background
		delete!(vis["Background"])
	end

    # Compute actor transformations
    actor_translations, actor_rotations = actor_transformations(solver, N)
	# We get rid of the last frame if the final state constraints are relaxed

	for k=1:N-1
        # Set the poses of the two actors.
        MeshCat.atframe(anim, vis, k) do frame
            for i=1:p
                settransform!(frame[scope*"/actor_bundle_$i"], compose(actor_translations[k][i],actor_rotations[k][i]))
            end
			# zoom = 1.0
			# setprop!(frame["/Cameras/default/rotated/<object>"], "zoom", zoom)
			setprop!(frame["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
            camera_transformation = compose(
                birdseye_trans,
                birdseye_rot)
            settransform!(frame["/Cameras/default"], camera_transformation)
        end
	end
	# setprop!(framevis["/Cameras/default/rotated/<object>"], "zoom", 0.5)
	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
	return vis, anim
end

function scene_animation(solver_leader::TO.ALTROSolver{T},
	solver_follower::TO.ALTROSolver{T},
	vis, plx; no_background=false) where T
	n,m,N = size(solver_leader)
	# Animate the scene
	birdseye_trans = Translation(0.00, 0.00, -2.45)
    birdseye_rot = compose(
		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))

	# Plot Trajectory
    anim = MeshCat.Animation()
	default_framerate = 3
	anim = MeshCat.Animation(anim.clips, default_framerate)
	delete!(vis["/Grid"])
	delete!(vis["/Axes"])
	if no_background
		delete!(vis["Background"])
	end
	delete!(vis["/Background"]) ####################################

    # Compute actor transformations
	actor_translations_leader, actor_rotations_leader = actor_transformations(solver_leader, N)
	actor_translations_follower, actor_rotations_follower = actor_transformations(solver_follower, N)
	# We get rid of the last frame if the final state constraints are relaxed


	alpha = 3e-1
	for k=1:N-1
        # Set the poses of the two actors.
        MeshCat.atframe(anim, vis, k) do frame
            settransform!(frame[scope*"/actor_bundle_1"],
				compose(actor_translations_leader[k][1],actor_rotations_leader[k][1]))
			settransform!(frame[scope*"/actor_bundle_2"],
				compose(actor_translations_follower[k][1],actor_rotations_follower[k][1]))
			# zoom = 1.0
			# setprop!(frame["/Cameras/default/rotated/<object>"], "zoom", zoom)
			setprop!(frame["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
            camera_transformation = compose(
                birdseye_trans,
                birdseye_rot)
            settransform!(frame["/Cameras/default"], camera_transformation)
        end
	end
	# setprop!(framevis["/Cameras/default/rotated/<object>"], "zoom", 0.5)
	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
	return vis, anim
end

function actor_transformations(solver::TO.ALTROSolver{T}, plx, iter) where T
    n,m,N = size(solver)
	X = TO.states(solver)
    actor_translations = []
    actor_rotations = []
    for k=1:iter
        actor_translation = [Translation(X[k][plx[i]]..., 0) for i=1:p]
        if k==N
            # Deals with the last step, we assume constant heading
            actor_rotation = actor_rotations[end]
        else
            angles = [atan(X[k+1][plx[i][2]]-X[k][plx[i][2]],
				X[k+1][plx[i][1]]-X[k][plx[i][1]]) for i=1:p]
            actor_rotation = [LinearMap(AngleAxis(angles[i], 0, 0, 1.)) for i=1:p]
        end
        push!(actor_translations, actor_translation)
        push!(actor_rotations, actor_rotation)
    end
    return actor_translations, actor_rotations
end
