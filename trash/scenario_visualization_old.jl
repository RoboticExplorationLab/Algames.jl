export
	solver_scope,
	animation,
	build_actors,
	build_trajectory,
	scene_animation,
	actor_transformations


function solver_scope(solver)
	return Symbol(typeof(solver).name)
end

# function solver_scope(solver)
# 	return Symbol(string(typeof(solver).name)*"_1")
# end

function animation(solver::TO.AbstractSolver, scenario::Scenario{T};
	vis::Visualizer=Visualizer(), anim::MeshCat.Animation=MeshCat.Animation(),
	open_vis::Bool=true,
	display_actors::Bool=true,
	display_trajectory::Bool=false,
	display_open_loop_plan::Bool=false,
	display_sampled_plan::Bool=false,
	no_background::Bool=false,
	actor_tracking::Bool=false,
	tracking_ind::Int=1,
	horizon::Int=size(solver)[3],
	dt::T=5e-1,
	scale::T=10.,
	camera_offset::Bool=true,
	camera_mvt::Bool=true,
	α_fading::T=0.9) where T

	# vis.core.scope
	# open_vis ? open(vis) : nothing
	delete!(vis["/Grid"])
	delete!(vis["/Axes"])
	no_background ? delete!(vis["Background"]) : nothing

	scope = solver_scope(solver)

	build_scenario(vis, scenario, scale=scale)
	if display_actors
		build_actors(vis, scenario, scope, scale=scale, α_fading=α_fading)
	elseif haskey(vis.core.tree.children, "meshcat") &&
		haskey(vis.core.tree.children["meshcat"].children, string(scope))
		delete!(vis.core.tree.children["meshcat"].children[string(solver_scope(solver))], "actors")
	end
	if display_trajectory
		build_trajectory(vis, solver, dt=dt, scale=scale, α_fading=α_fading)
	elseif haskey(vis.core.tree.children, "meshcat") &&
		haskey(vis.core.tree.children["meshcat"].children, string(scope))
		delete!(vis.core.tree.children["meshcat"].children[string(scope)], "trajectory")
	end

	if display_open_loop_plan && solver_scope(solver) == :UKFGamesSolver
		build_open_loop_plan(vis, solver, horizon=horizon, dt=dt, scale=scale, α_fading=α_fading)
	elseif haskey(vis.core.tree.children, "meshcat") &&
		haskey(vis.core.tree.children["meshcat"].children, string(scope))
		delete!(vis.core.tree.children["meshcat"].children[string(scope)], "open_loop_plan")
	end

	if display_sampled_plan && solver_scope(solver) == :UKFGamesSolver
		build_sampled_plan(vis, solver, horizon=horizon, dt=dt, scale=scale, α_fading=α_fading)
	elseif haskey(vis.core.tree.children, "meshcat") &&
		haskey(vis.core.tree.children["meshcat"].children, string(scope))
		delete!(vis.core.tree.children["meshcat"].children[string(scope)], "sampled_plan")
	end

	# Animate the scene
	default_framerate = 3
	anim = MeshCat.Animation(anim.clips, default_framerate)
	scene_animation(solver, vis, anim,
		actor_tracking=actor_tracking,
		tracking_ind=tracking_ind,
		horizon=horizon,
		dt=dt,
		scale=scale,
		camera_offset=camera_offset,
		camera_mvt=camera_mvt,
		)
    MeshCat.setanimation!(vis, anim)
	open_vis ? open(vis) : nothing
	return vis, anim
end

function build_actors(vis::Visualizer, scenario::Scenario{T}, scope::Symbol; scale::T=10., α_fading::T=0.9) where T
	actors_radii = scenario.actors_radii*scale
	actors_types = scenario.actors_types
	# car_offset = compose(compose(compose(
	# 		Translation(0.0, 0.0, -0.0),
	# 		LinearMap(AngleAxis(pi/200, 0, 1, 0))),
	# 		LinearMap(AngleAxis(pi/2, 0, 0, 1))),
	# 		LinearMap(AngleAxis(pi/2, 1, 0, 0)))

	car_offset = compose(compose(compose(
			Translation(-0.03, 0.0, 0.31),
			LinearMap(AngleAxis(-0.012+pi/200, 0, 1, 0))),
			LinearMap(AngleAxis(pi/2, 0, 0, 1))),
			LinearMap(AngleAxis(pi/2, 1, 0, 0)))
	ped_offset = LinearMap(AngleAxis(pi/2, 0, 0, 1))

	cylinder_height_car = 0.01*scale
	cylinder_height_ped = 0.035*scale
	color_values = [100/255 149/255 237/255; # cornflowerblue
					# 100/255 149/255 237/255;       # orange
					# 100/255 149/255 237/255;       # orange
					1       125/255 0;       # orange
					1       125/255 0;       # orange
					100/255 149/255 237/255;       # orange
					34/255  139/255 34/255;  # forestgreen
					1       0       0;       # red
					1       1       0;       # yellow
					0.5     0.5     0.5;     # gray
					1       125/255 0;       # orange
					0       1       0]       # Lime
	colors = ["cornflowerblue",
		# "cornflowerblue",
		# "cornflowerblue",
		"orange",
		"orange",
		"cornflowerblue",
		"forestgreen", "red", "yellow",
		 "gray", "orange", "lime"]

	pkg_path = joinpath(dirname(@__FILE__), "../../")

	path = String(scope)*"/actors/actor_bundle_"
	for i = 1:length(actors_radii)
		if actors_types[i] == :car
			# actor = load(joinpath(pkg_path, "resources/object/car_geometry.obj"), GLUVMesh)
			# actor.vertices .*= 0.025*scale
			# actor_image = PngImage(joinpath(pkg_path, "resources/textures/"*colors[i]*".png"))
			# actor_texture = Texture(image=actor_image)
			# actor_material = MeshLambertMaterial(map=actor_texture)
			# setobject!(vis[path*"$i/actor"], actor, actor_material)
			obj_path = joinpath(pkg_path, "resources/object/car_geometry.obj")
			mtl_path = joinpath(pkg_path, "resources/material/"*colors[i]*".mtl")
			actor = ModifiedMeshFileObject(obj_path, mtl_path, scale=0.106)
			setobject!(vis[path*"$i/actor"], actor)
			settransform!(vis[path*"$i/actor"], car_offset)

			# Add collision avoidance cylinders
			collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,cylinder_height_car), actors_radii[i]) ### need to plaot cylinders per timestep
			cylinder_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., α_fading))
			setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)

		elseif actors_types[i] == :pedestrian
			actor = load(joinpath(pkg_path, "resources/object/pedestrian_geometry.obj"), GLUVMesh)
			actor.vertices .*= 0.001/15 * scale
			# actor_image = PngImage(joinpath(pkg_path, "resources/textures/light_boundary.png"))
			# actor_texture = Texture(image=actor_image)
			# actor_material = MeshLambertMaterial(map=actor_texture)
			actor_material = MeshPhongMaterial(color=RGBA(0., 0., 0., 1.0))
			setobject!(vis[path*"$i/actor"], actor, actor_material)
			settransform!(vis[path*"$i/actor"], ped_offset)

			# Add collision avoidance cylinders
			collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,cylinder_height_ped), actors_radii[i]) ### need to plaot cylinders per timestep
			cylinder_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., α_fading))
			setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)
		end
	end
	return nothing
end

function build_trajectory(vis::Visualizer, solver::TO.AbstractSolver{T}; dt::T=5e-1, scale::T=10., α_fading::T=0.9) where {T}
	build_trajectory(vis,solver,TO.states(solver),dt=dt,scale=scale,α_fading=α_fading)
end

function build_trajectory(vis::Visualizer, solver::TO.AbstractSolver{T}, X::Vector{K}; dt::T=5e-1, scale::T=10., α_fading::T=0.9) where {T,K}
	color_values = [100/255 149/255 237/255; # cornflowerblue
					# 100/255 149/255 237/255;       # orange
					# 100/255 149/255 237/255;       # orange
					1       125/255 0;       # orange
					1       125/255 0;       # orange
					# 100/255 149/255 237/255;       # orange
					34/255  139/255 34/255;  # forestgreen
					1       0       0;       # red
					1       1       0;       # yellow
					0.5     0.5     0.5;     # gray
					1       125/255 0;       # orange
					0       1       0]       # Lime
	colors = ["cornflowerblue",
		# "cornflowerblue",
		# "cornflowerblue",
		"orange",
		"orange",
		# "cornflowerblue",
		"forestgreen", "red", "yellow",
		 "gray", "orange", "lime"]
	cyl_height = 0.005*scale
	cyl_radius = 0.015*scale

	n,m,N = size(solver)
	model = TO.get_model(solver)
	p = model.p
	px = model.px
	path = String(solver_scope(solver))*"/trajectory"
	for k = 1:length(X)
		for i = 1:p
			cyl_pos = Translation(scale*X[k][px[i]]..., 0)
			aug_cyl_height = cyl_height + i*0.00 + k*0.01/length(X)*scale
			cyl = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,aug_cyl_height), cyl_radius)
			cyl_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., α_fading))
			setobject!(vis[path*"/cyl_$i/step"*"_$k"], cyl, cyl_material)
			settransform!(vis[path*"/cyl_$i/step"*"_$k"], cyl_pos)
		end
	end
	return nothing
end

function build_open_loop_plan(vis::Visualizer, solver::TO.AbstractSolver{T};
	horizon::Int=solver.solver[1].N, dt::T=5e-1, scale::T=10., α_fading::T=0.9) where {T,K}

	color_values = [100/255 149/255 237/255; # cornflowerblue
					100/255 149/255 237/255; # cornflowerblue
					1       125/255 0;       # orange
					34/255  139/255 34/255;  # forestgreen
					1       0       0;       # red
					1       1       0;       # yellow
					0.5     0.5     0.5;     # gray
					1       125/255 0;       # orange
					0       1       0]       # Lime
	colors = ["cornflowerblue",
		"cornflowerblue",
		"orange", "forestgreen", "red", "yellow",
		"gray", "orange", "lime"]
	cyl_height = 0.007*scale
	cyl_radius = 0.010*scale

	n,m,N = size(solver)
	model = TO.get_model(solver)
	n,m,pu,p = size(model)
	px = model.px
	path = String(solver_scope(solver))*"/open_loop_plan"

	for j = 1:horizon
		for i = 1:p
			aug_cyl_height = cyl_height# + i*0.00 + k*0.01/M*scale
			cyl = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,aug_cyl_height), cyl_radius)
			cyl_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., α_fading))
			setobject!(vis[path*"/cyl_$i/stage_$j"], cyl, cyl_material)
		end
	end
	return nothing
end

function build_sampled_plan(vis::Visualizer, solver::TO.AbstractSolver{T};
	horizon::Int=solver.solver[1].N, dt::T=5e-1, scale::T=10., α_fading::T=0.9) where {T,K}

	color_values = [100/255 149/255 237/255; # cornflowerblue
					1       125/255 0;       # orange
					34/255  139/255 34/255;  # forestgreen
					1       0       0;       # red
					1       1       0;       # yellow
					0.5     0.5     0.5;     # gray
					1       125/255 0;       # orange
					0       1       0]       # Lime
	colors = ["cornflowerblue", "orange", "forestgreen", "red", "yellow",
		 "gray", "orange", "lime"]
	cyl_height = 0.020*scale
	cyl_radius = 0.005*scale

	n,m,N = size(solver)
	model = TO.get_model(solver)
	n,m,pu,p = size(model)
	px = model.px
	path = String(solver_scope(solver))*"/sampled_plan"
	for l = 1:solver.r-2
		for j = 1:horizon
			for i = 1:p
				aug_cyl_height = cyl_height# + i*0.00 + k*0.01/M*scale
				cyl = Cylinder(Point(0.0,0.0,0.0),
					Point(0.0,0.0,aug_cyl_height), cyl_radius)
				cyl_material = MeshPhongMaterial(color=RGBA(color_values[4,:]..., α_fading))
				setobject!(vis[path*"/cyl_$i/stage_$j/sample_$l"], cyl, cyl_material)
			end
		end
	end
	return nothing
end

function scene_animation(solver::TO.AbstractSolver, vis::Visualizer,
	anim::MeshCat.Animation; actor_tracking::Bool=false, tracking_ind::Int=1,
	horizon::Int=solver.solver[1].N, dt::T=5e-1, scale::T=10., camera_offset::Bool=true, camera_mvt::Bool=true) where T
	n,m,N = size(solver)
	model = TO.get_model(solver)
	p = model.p
	px = model.px
	offset = camera_offset ? -0.80*scale : 0.00
	###
	# offset = -0.00*scale
	###
	# Animate the scene
	birdseye_trans = Translation(0.0, 0.0, 1.9*scale)
    birdseye_rot = compose(
		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))
	###
	# birdseye_trans = Translation(0.0, 0.0, 1.0*scale)
    # birdseye_rot = compose(
	# 	LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
	# 	LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))
	###

    # Compute actor transformations
    actor_translations, actor_rotations = actor_transformations(solver, dt=dt, scale=scale)
	if solver_scope(solver) == :UKFGamesSolver
		start_ind, end_ind, α = time_resampling([z.dt for z in solver.stats.traj], dt=dt)
	end

	# We get rid of the last frame if the final state constraints are relaxed
	actor_path = String(solver_scope(solver))*"/actors/actor_bundle"
	trajectory_path = String(solver_scope(solver))*"/trajectory"
	open_loop_plan_path = String(solver_scope(solver))*"/open_loop_plan"
	sampled_plan_path = String(solver_scope(solver))*"/sampled_plan"
	roadway_path = "roadway"
	for k=1:length(actor_translations)-1
        # Set the poses of the two actors.
        MeshCat.atframe(anim, k) do
			mean_x_dot = camera_mvt ? mean([actor_translations[k][j].translation[1] for j=1:p]) : 0.0
			roadway_translation = Translation([offset-mean_x_dot, 0., 0.])
			settransform!(
					vis[roadway_path],
					roadway_translation)
			settransform!(
					vis[trajectory_path],
					roadway_translation)

			if solver_scope(solver) == :UKFGamesSolver
				for j=1:horizon
					x_start = TO.state(solver.stats.traj_pred[start_ind[k]][j])
					x_end = TO.state(solver.stats.traj_pred[end_ind[k]][j])
					x = x_start + α[k]*(x_end-x_start)
					for i=1:p
						cyl_pos = compose(roadway_translation, Translation(scale*x[px[i]]..., 0))
						settransform!(vis[open_loop_plan_path*"/cyl_$i/stage_$j"], cyl_pos)
					end
				end
				if k == 1
					for j=1:horizon
						x_start = TO.state(solver.stats.traj_pred[start_ind[k]][j])
						x_end = TO.state(solver.stats.traj_pred[end_ind[k]][j])
						x = x_start + α[k]*(x_end-x_start)
						for i=1:p
							cyl_pos = compose(roadway_translation, Translation(scale*x[px[i]]..., 0))
							for l = 1:solver.r-2
								settransform!(vis[sampled_plan_path*"/cyl_$i/stage_$j/sample_$l"], cyl_pos)
							end
						end
					end
				else
					for l = 1:solver.r-2
						for j=1:horizon
							x_start = TO.state(solver.stats.traj_sampled[start_ind[k]+1][l][j])
							x_end = TO.state(solver.stats.traj_sampled[end_ind[k]+1][l][j])
							x = x_start + α[k]*(x_end-x_start) #+ rand(length(x_start))

							for i=1:p
								cyl_pos = compose(roadway_translation, Translation(scale*x[px[i]]..., 0))
								settransform!(vis[sampled_plan_path*"/cyl_$i/stage_$j/sample_$l"], cyl_pos)
							end
						end
					end
				end
			end

            for i=1:p
                settransform!(
					vis[actor_path*"_$i"],
					compose(compose(
						roadway_translation,
						actor_translations[k][i]),
						actor_rotations[k][i]))
            end
			setprop!(vis["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
			setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 1.0)
			camera_transformation = compose(
				birdseye_trans,
				birdseye_rot)
            # camera_transformation = compose(compose(
			# 	Translation(-4.0, 0.0, -k*0.6),
			# 	compose(
			# 	birdseye_trans,
			# 	birdseye_rot)),Translation(0.0, -k*0.3, -5+k*0.6))
			settransform!(vis["/Cameras/default"], camera_transformation)

			# setprop!(
			# 	# vis["/Cameras/default"],
			# 	vis["/Cameras/default/rotated/<object>"],
			# 	"position",
			# 	[-15.0 + 1.0*k, 0., 00.])
        end
	end
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 3.0)
	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
	return nothing
end


function actor_transformations(solver::TO.AbstractSolver{T}; dt::T=5e-1, scale::T=10.) where {T}
	actor_transformations(solver,TO.states(solver),dt=dt,scale=scale)
end

function actor_transformations(solver::TO.AbstractSolver{T}, X::Vector{K}; dt::T=5e-1, scale::T=10.) where {T,K}
	n,m,N = size(solver)
	px = TO.get_model(solver).px
	p = length(px)
	# X = TO.states(solver)
    actor_translations = []
    actor_rotations = []
	iter = length(X)
    for k=1:iter-1
        actor_translation = [Translation(scale*X[k][px[i]]..., 0) for i=1:p]
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

function state_resampling(traj::Vector{K}; dt::T=5e-1) where {T,K}
	n = length(traj[1]._x)
	start_ind, end_ind, α = time_resampling([z.dt for z in traj], dt=dt)
	M = length(α)
	# Δt = sum([z.dt for z in traj])
	# sum_dt = cumsum([0.;[z.dt for z in traj]])
	# M = Int(floor(Δt/dt)-1)
	X = Vector{SVector{n,T}}(undef, M)
	# t = 0.
	for k = 1:M
		# start_ind = findlast(x -> x <= t, sum_dt)
		# end_ind = start_ind+1
		x_start = TO.state(traj[start_ind[k]])
		x_end = TO.state(traj[end_ind[k]])
		# α = (t - sum_dt[start_ind])/(sum_dt[end_ind] - sum_dt[start_ind])
		X[k] = x_start + α[k]*(x_end-x_start)
		# t += dt
	end
	return X
end

function time_resampling(dts::Vector{T}; dt::T=5e-1) where {T}
	Δt = sum(dts)
	sum_dt = cumsum([0.;dts])
	M = Int(floor(Δt/dt)-1)
	t = 0.
	start_ind = zeros(Int, M)
	end_ind = zeros(Int, M)
	α = zeros(T, M)
	for k = 1:M
		start_ind[k] = findlast(x -> x <= t, sum_dt)
		end_ind[k] = start_ind[k]+1
		α[k] = (t - sum_dt[start_ind[k]])/(sum_dt[end_ind[k]] - sum_dt[start_ind[k]])
		t += dt
	end
	return start_ind, end_ind, α
end



function ModifiedMeshFileObject(obj_path::String, material_path::String; scale::T=0.1) where {T}
    obj = MeshFileObject(obj_path)
    rescaled_contents = rescale_contents(obj_path, scale=scale)
    material = select_material(material_path)
    mod_obj = MeshFileObject(
        rescaled_contents,
        obj.format,
        material,
        obj.resources,
        )
    return mod_obj
end

function rescale_contents(obj_path::String; scale::T=0.1) where T
    lines = readlines(obj_path)
    rescaled_lines = copy(lines)
    for (k,line) in enumerate(lines)
        if length(line) >= 2
            if line[1] == 'v'
                stringvec = split(line, " ")
                vals = map(x->parse(Float64,x),stringvec[end-2:end])
                rescaled_vals = vals .* scale
                rescaled_lines[k] = join([stringvec[1]; string.(rescaled_vals)], " ")
            end
        end
    end
    rescaled_contents = join(rescaled_lines, "\r\n")
    return rescaled_contents
end

function select_material(material_path::String)
    mtl_file = open(material_path)
    mtl = read(mtl_file, String)
    return mtl
end










# function animation(solver::TO.AbstractSolver, scenario::Scenario{T};
# 	vis::Visualizer=Visualizer(), anim::MeshCat.Animation=MeshCat.Animation(),
# 	open_vis::Bool=true,
# 	display_actors::Bool=true,
# 	display_trajectory::Bool=false,
# 	no_background::Bool=false,
# 	actor_tracking::Bool=false,
# 	tracking_ind::Int=1,
# 	scale::T=10.) where T
#
# 	# vis.core.scope
# 	# open_vis ? open(vis) : nothing
# 	delete!(vis["/Grid"])
# 	delete!(vis["/Axes"])
# 	no_background ? delete!(vis["Background"]) : nothing
#
# 	scope = solver_scope(solver)
#
# 	build_scenario(vis, scenario, scale=scale)
# 	if display_actors
# 		build_actors(vis, scenario, scope, scale=scale)
# 	elseif haskey(vis.core.tree.children, "meshcat") &&
# 		haskey(vis.core.tree.children["meshcat"].children, string(scope))
# 		delete!(vis.core.tree.children["meshcat"].children[string(solver_scope(solver))], "actors")
# 	end
# 	if display_trajectory
# 		build_trajectory(vis, solver)
# 	elseif haskey(vis.core.tree.children, "meshcat") &&
# 		haskey(vis.core.tree.children["meshcat"].children, string(scope))
# 		delete!(vis.core.tree.children["meshcat"].children[string(scope)], "trajectory")
# 	end
#
# 	# Animate the scene
# 	default_framerate = 3
# 	anim = MeshCat.Animation(anim.clips, default_framerate)
# 	scene_animation(solver, vis, anim,
# 		actor_tracking=actor_tracking,
# 		tracking_ind=tracking_ind,
# 		scale=scale)
#     MeshCat.setanimation!(vis, anim)
# 	open_vis ? open(vis) : nothing
# 	return vis, anim
# end
#
# function build_actors(vis::Visualizer, scenario::Scenario{T}, scope::Symbol; scale::T=10.) where T
# 	actors_radii = scenario.actors_radii*scale
# 	actors_types = scenario.actors_types
# 	car_offset = compose(compose(compose(
# 			Translation(0.0, 0.0, -0.0),
# 			LinearMap(AngleAxis(pi/200, 0, 1, 0))),
# 			LinearMap(AngleAxis(pi/2, 0, 0, 1))),
# 			LinearMap(AngleAxis(pi/2, 1, 0, 0)))
# 	ped_offset = LinearMap(AngleAxis(pi/2, 0, 0, 1))
#
# 	cylinder_height_car = 0.01*scale
# 	cylinder_height_ped = 0.035*scale
# 	color_values = [100/255 149/255 237/255; # cornflowerblue
# 					1       125/255 0;       # orange
# 					34/255  139/255 34/255;  # forestgreen
# 					1       0       0;       # red
# 					1       1       0;       # yellow
# 					0.5     0.5     0.5;     # gray
# 					1       165/255 0;       # orange
# 					0       1       0]       # Lime
# 	colors = ["cornflowerblue", "orange", "forestgreen", "red", "yellow",
# 		 "gray", "orange", "lime"]
#
# 	pkg_path = joinpath(dirname(@__FILE__), "../../")
#
# 	path = String(scope)*"/actors/actor_bundle_"
# 	for i = 1:length(actors_radii)
# 		if actors_types[i] == :car
# 			actor = load(joinpath(pkg_path, "resources/object/car_geometry.obj"), GLUVMesh)
# 			actor.vertices .*= 0.025*scale
# 			actor_image = PngImage(joinpath(pkg_path, "resources/textures/"*colors[i]*".png"))
# 			actor_texture = Texture(image=actor_image)
# 			actor_material = MeshLambertMaterial(map=actor_texture)
# 			setobject!(vis[path*"$i/actor"], actor, actor_material)
# 			settransform!(vis[path*"$i/actor"], car_offset)
#
# 			# Add collision avoidance cylinders
# 			collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
# 				Point(0.0,0.0,cylinder_height_car), actors_radii[i]) ### need to plaot cylinders per timestep
# 			cylinder_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., 9e-1))
# 			setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)
#
# 		elseif actors_types[i] == :pedestrian
# 			actor = load(joinpath(pkg_path, "resources/object/pedestrian_geometry.obj"), GLUVMesh)
# 			actor.vertices .*= 0.001/15 * scale
# 			actor_image = PngImage(joinpath(pkg_path, "resources/textures/black_boundary.png"))
# 			actor_texture = Texture(image=actor_image)
# 			actor_material = MeshLambertMaterial(map=actor_texture)
# 			setobject!(vis[path*"$i/actor"], actor, actor_material)
# 			settransform!(vis[path*"$i/actor"], ped_offset)
#
# 			# Add collision avoidance cylinders
# 			collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
# 				Point(0.0,0.0,cylinder_height_ped), actors_radii[i]) ### need to plaot cylinders per timestep
# 			cylinder_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., 6e-1))
# 			setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)
# 		end
# 	end
# 	return nothing
# end
#
# function build_trajectory(vis::Visualizer, solver::TO.AbstractSolver{T}; scale::T=10.) where {T}
# 	color_values = [100/255 149/255 237/255; # cornflowerblue
# 					1       125/255 0;       # orange
# 					34/255  139/255 34/255;  # forestgreen
# 					1       0       0;       # red
# 					1       1       0;       # yellow
# 					0.5     0.5     0.5;     # gray
# 					1       165/255 0;       # orange
# 					0       1       0]       # Lime
# 	colors = ["cornflowerblue", "orange", "forestgreen", "red", "yellow",
# 		 "gray", "orange", "lime"]
# 	cyl_height = 0.005*scale
# 	cyl_radius = 0.025*scale
#
# 	n,m,N = size(solver)
# 	model = TO.get_model(solver)
# 	n,m,pu,p = size(model)
# 	px = model.px
# 	path = String(solver_scope(solver))*"/trajectory"
# 	X = TO.states(solver)
# 	for k = 1:length(X)
# 		for i = 1:p
# 			cyl_pos = Translation(scale*X[k][px[i]]..., 0)
# 			aug_cyl_height = cyl_height + i*0.00 + k*0.01/length(X)*scale
# 			cyl = Cylinder(Point(0.0,0.0,0.0),
# 				Point(0.0,0.0,aug_cyl_height), cyl_radius)
# 			cyl_material = MeshPhongMaterial(color=RGBA(color_values[i,:]..., 9e-1))
# 			setobject!(vis[path*"/cyl_$i/step"*"_$k"], cyl, cyl_material)
# 			settransform!(vis[path*"/cyl_$i/step"*"_$k"], cyl_pos)
# 		end
# 	end
# 	return nothing
# end
#
# function scene_animation(solver::TO.AbstractSolver, vis::Visualizer,
# 	anim::MeshCat.Animation; actor_tracking::Bool=false, tracking_ind::Int=1, scale::T=10.) where T
# 	n,m,N = size(solver)
# 	model = TO.get_model(solver)
# 	n,m,pu,p = size(model)
# 	X = TO.states(solver)
# 	U = TO.controls(solver)
#
# 	len = length(solver.Z)
#
# 	# Animate the scene
# 	birdseye_trans = Translation(0.0, 0.0, 0.0*scale)
# 	# birdseye_trans = Translation(0.0, 0.0, 0.5*scale)
# 	# birdseye_trans = Translation(0.0, 0.0, 1.1*scale)
# 	# birdseye_trans = Translation(0.0, 0.0, 2.2)
#     birdseye_rot = compose(
# 		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
# 		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))
# 		# LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
# 		# LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))
#     # Compute actor transformations
#     actor_translations, actor_rotations = actor_transformations(solver, len, scale=scale)
# 	# We get rid of the last frame if the final state constraints are relaxed
# 	actor_path = String(solver_scope(solver))*"/actors/actor_bundle"
# 	roadway_path = "roadway"
# 	trajectory_path = String(solver_scope(solver))*"/trajectory"
# 	for k=1:len-1
#         # Set the poses of the two actors.
#         MeshCat.atframe(anim, k) do
# 			mean_x_dot = mean([actor_translations[k][j].translation[1] for j=1:p])
# 			roadway_translation = Translation([-mean_x_dot, 0., 0.])
#             for i=1:p
#                 settransform!(
# 					vis[actor_path*"_$i"],
# 					compose(compose(
# 						roadway_translation,
# 						actor_translations[k][i]),
# 						actor_rotations[k][i]))
# 				settransform!(
# 						vis[roadway_path],
# 						roadway_translation)
# 				settransform!(
# 						vis[trajectory_path],
# 						roadway_translation)
#             end
# 			setprop!(vis["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
# 			setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 1.0)
#             camera_transformation = compose(
# 				birdseye_trans,
#                 birdseye_rot)
# 			settransform!(vis["/Cameras/default"], camera_transformation)
#
# 			# setprop!(
# 			# 	# vis["/Cameras/default"],
# 			# 	vis["/Cameras/default/rotated/<object>"],
# 			# 	"position",
# 			# 	[-15.0 + 1.0*k, 0., 00.])
#         end
# 	end
# 	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", max(1.0,5.0/scale))
# 	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
# 	return nothing
# end
#
# function actor_transformations(solver, iter; scale::T=10.) where T
# 	n,m,N = size(solver)
# 	px = TO.get_model(solver).px
# 	p = length(px)
# 	X = TO.states(solver)
#     actor_translations = []
#     actor_rotations = []
#     for k=1:iter
#         actor_translation = [Translation(scale*X[k][px[i]]..., 0) for i=1:p]
#         if k==iter
#             # Deals with the last step, we assume constant heading
#             actor_rotation = actor_rotations[end]
#         else
#             angles = [atan(X[k+1][px[i][2]]-X[k][px[i][2]],
# 					  X[k+1][px[i][1]]-X[k][px[i][1]]) for i=1:p]
#             actor_rotation = [LinearMap(AngleAxis(angles[i], 0, 0, 1.)) for i=1:p]
#         end
#         push!(actor_translations, actor_translation)
#         push!(actor_rotations, actor_rotation)
#     end
#     return actor_translations, actor_rotations
# end




################################
################################
################################
#
#
# function animation(solver_leader::TO.ALTROSolver{T},
# 	solver_follower::TO.ALTROSolver{T},
# 	scenario::Scenario{T}, plx;
# 	no_background=false) where T
#
# 	n,m,N = size(solver_leader)
#     # Open visualizer
#     vis = Visualizer()
# 	vis.core.scope
#     open(vis)
# 	build_actors(vis, scenario, scale=scale)
# 	build_scenario(vis, scenario, scale=scale)
# 	# Animate the scene
# 	vis, anim = scene_animation(solver_leader, solver_follower,
# 		vis, plx)
#     MeshCat.setanimation!(vis, anim)
# end
#
# function animation(solver::TO.ALTROSolver{T}, scenario::Scenario{T}, plx;
# 	no_background=false) where T
#
# 	n,m,N = size(solver)
#     # Open visualizer
#     vis = Visualizer()
# 	vis.core.scope
#     open(vis)
# 	build_actors(vis, scenario, scale=scale)
# 	build_scenario(vis, scenario, scale=scale)
# 	# Animate the scene
# 	vis, anim = scene_animation(solver, vis, plx; no_background=no_background)
#     MeshCat.setanimation!(vis, anim)
# end
#
# function scene_animation(solver::TO.ALTROSolver{T}, vis, plx; no_background=false) where T
# 	n,m,N = size(solver)
# 	X = TO.states(solver)
# 	U = TO.controls(solver)
#
# 	# Animate the scene
# 	birdseye_trans = Translation(0.0, 0.0, -2.2)
#     birdseye_rot = compose(
# 		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
# 		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))
#
# 	# Plot Trajectory
#     anim = MeshCat.Animation()
# 	default_framerate = 3
# 	anim = MeshCat.Animation(anim.clips, default_framerate)
# 	delete!(vis["/Grid"])
# 	delete!(vis["/Axes"])
# 	if no_background
# 		delete!(vis["Background"])
# 	end
#
#     # Compute actor transformations
#     actor_translations, actor_rotations = actor_transformations(solver, N, scale=scale)
# 	# We get rid of the last frame if the final state constraints are relaxed
#
# 	for k=1:N-1
#         # Set the poses of the two actors.
#         MeshCat.atframe(anim, vis, k) do frame
#             for i=1:p
#                 settransform!(frame[scope*"/actor_bundle_$i"], compose(actor_translations[k][i],actor_rotations[k][i]))
#             end
# 			# zoom = 1.0
# 			# setprop!(frame["/Cameras/default/rotated/<object>"], "zoom", zoom)
# 			setprop!(frame["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
#             camera_transformation = compose(
#                 birdseye_trans,
#                 birdseye_rot)
#             settransform!(frame["/Cameras/default"], camera_transformation)
#         end
# 	end
# 	# setprop!(framevis["/Cameras/default/rotated/<object>"], "zoom", 0.5)
# 	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
# 	return vis, anim
# end
#
# function scene_animation(solver_leader::TO.ALTROSolver{T},
# 	solver_follower::TO.ALTROSolver{T},
# 	vis, plx; no_background=false) where T
# 	n,m,N = size(solver_leader)
# 	# Animate the scene
# 	birdseye_trans = Translation(0.00, 0.00, -2.45)
#     birdseye_rot = compose(
# 		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
# 		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))
#
# 	# Plot Trajectory
#     anim = MeshCat.Animation()
# 	default_framerate = 3
# 	anim = MeshCat.Animation(anim.clips, default_framerate)
# 	delete!(vis["/Grid"])
# 	delete!(vis["/Axes"])
# 	if no_background
# 		delete!(vis["Background"])
# 	end
# 	delete!(vis["/Background"]) ####################################
#
#     # Compute actor transformations
# 	actor_translations_leader, actor_rotations_leader = actor_transformations(solver_leader, N, scale=scale)
# 	actor_translations_follower, actor_rotations_follower = actor_transformations(solver_follower, N, scale=scale)
# 	# We get rid of the last frame if the final state constraints are relaxed
#
#
# 	alpha = 3e-1
# 	for k=1:N-1
#         # Set the poses of the two actors.
#         MeshCat.atframe(anim, vis, k) do frame
#             settransform!(frame[scope*"/actor_bundle_1"],
# 				compose(actor_translations_leader[k][1],actor_rotations_leader[k][1]))
# 			settransform!(frame[scope*"/actor_bundle_2"],
# 				compose(actor_translations_follower[k][1],actor_rotations_follower[k][1]))
# 			# zoom = 1.0
# 			# setprop!(frame["/Cameras/default/rotated/<object>"], "zoom", zoom)
# 			setprop!(frame["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
#             camera_transformation = compose(
#                 birdseye_trans,
#                 birdseye_rot)
#             settransform!(frame["/Cameras/default"], camera_transformation)
#         end
# 	end
# 	# setprop!(framevis["/Cameras/default/rotated/<object>"], "zoom", 0.5)
# 	# setprop!(vis, "/Lights/DirectionalLight/<object>", "intensity", 1.2)
# 	return vis, anim
# end
#
# function actor_transformations(solver::TO.ALTROSolver{T}, plx, iter) where T
#     n,m,N = size(solver)
# 	X = TO.states(solver)
#     actor_translations = []
#     actor_rotations = []
#     for k=1:iter
#         actor_translation = [Translation(X[k][plx[i]]..., 0) for i=1:p]
#         if k==N
#             # Deals with the last step, we assume constant heading
#             actor_rotation = actor_rotations[end]
#         else
#             angles = [atan(X[k+1][plx[i][2]]-X[k][plx[i][2]],
# 				X[k+1][plx[i][1]]-X[k][plx[i][1]]) for i=1:p]
#             actor_rotation = [LinearMap(AngleAxis(angles[i], 0, 0, 1.)) for i=1:p]
#         end
#         push!(actor_translations, actor_translation)
#         push!(actor_rotations, actor_rotation)
#     end
#     return actor_translations, actor_rotations
# end
