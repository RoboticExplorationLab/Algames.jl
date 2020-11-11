export
	solver_scope,
	animation,
	build_actors,
	build_trajectory,
	scene_animation,
	actor_transformations,
	AnimationOptions,
	delete_object!,
	animation_colors


mutable struct Ellipsoid{n,T}
	mean::SVector{n,T}
	cov::Matrix{T}
end

@with_kw mutable struct AnimationOptions{T}
    # Options
	"Suffix added at the end of the scope of the solver."
	scope_suffix::String=""

	"Display the actors along with their collision avoidance cylinders."
	display_actors::Bool=true

	"Display the ground-truth trajectories of each actor."
	display_trajectory::Bool=false

	"Display the open-loop plan computed by the human-driven vehicle for all the actors."
	display_open_loop_plan::Bool=false

	"Display the open-loop plan sampled by the human-driven vehicle for all the actors."
	display_sampled_plan::Bool=false

	"Display ellipse contour corresponding to the uncertainty about the position of human driver vehicles."
	display_ellipse_contour::Bool=false

	"Number of points on the border of the ellipse contour."
	num_contour_points::Int=100

	"Hide the blue background."
	no_background::Bool=false

	"The camera tracks one actor."
	actor_tracking::Bool=false

	"The ID of the actor tracked by the camera."
	tracking_ind::Int=1

	"Number of time steps of the open-loop plans."
	horizon::Int=10#size(solver)[3]

	"Resampling time step"
	dt::T=25e-2

	"Geometric rescaling of the scene"
	scale::T=10.

	"Offset the camera to better see the open-loop plans."
	camera_offset::Bool=true

	"Allow the camera to follow the group of actors."
	camera_mvt::Bool=true

	"Fading to the trajectories and the collision cylinders."
	α_fading::T=0.9

	"Framerate of the animation."
	framerate::Int=3

	"Set of colors used to mark each actor."
	color_set::Vector{String}=["orange", "cornflowerblue", "forestgreen", "red", "yellow"]
end

function AnimationOptions(solver::TO.AbstractSolver{T}) where T
	n,m,N = size(solver)
	opts = AnimationOptions{T}(horizon=N)
	return opts
end

function solver_scope(solver, opts::AnimationOptions)
	str = string(typeof(solver).name) * opts.scope_suffix
	return Symbol(str)
end

function delete_object!(vis::Visualizer, scope::Symbol, object::String)
	if haskey(vis.core.tree.children, "meshcat") &&
		haskey(vis.core.tree.children["meshcat"].children, string(scope))
		delete!(vis.core.tree.children["meshcat"].children[string(scope)], object)
	end
	return nothing
end

function animation_colors(col::String)
	color_names = [
		"cornflowerblue",
		"orange",
		"forestgreen",
		"red",
		"yellow",
		"gray",
		"orange",
		"lime"]

	color_values = [
		100/255 149/255 237/255; # cornflowerblue
		1       125/255 0;       # orange
		34/255  139/255 34/255;  # forestgreen
		1       0       0;       # red
		1       1       0;       # yellow
		0.5     0.5     0.5;     # gray
		1       125/255 0;       # orange
		0       1       0]       # lime
	ind = findfirst(x -> x == col, color_names)
	if ind == nothing
		@show "Unknown color name."
	end
	return color_values[ind,:]
end

function animation_colors(index::Int, opts::AnimationOptions{T}) where {T}
	color_name = opts.color_set[index]
	return animation_colors(color_name)
end

function animation(solver::TO.AbstractSolver, scenario::Scenario{T};
	vis::Visualizer=Visualizer(), anim::MeshCat.Animation=MeshCat.Animation(),
	opts::AnimationOptions{T}=AnimationOptions(solver)) where T

	scope = solver_scope(solver, opts)
	delete!(vis["/Grid"])
	delete!(vis["/Axes"])
	opts.no_background ? delete!(vis["Background"]) : nothing

	build_scenario(vis, scenario, scale=opts.scale)
	if opts.display_actors
		build_actors(vis, scenario, scope, opts)
	else
		delete_object!(vis,scope,"actors")
	end
	if opts.display_trajectory
		build_trajectory(vis, solver, opts)
	else
		delete_object!(vis,scope,"trajectory")
	end

	if opts.display_open_loop_plan && occursin("UKFGamesSolver", string(scope))
		build_open_loop_plan(vis, solver, opts)
	else
		delete_object!(vis,scope,"open_loop_plan")
	end

	if opts.display_sampled_plan && occursin("UKFGamesSolver", string(scope))
		build_sampled_plan(vis, solver, opts)
	else
		delete_object!(vis,scope,"sampled_plan")
	end

	if opts.display_ellipse_contour && occursin("UKFGamesSolver", string(scope))
		build_ellipse_contour(vis, solver, scenario, opts)
	else
		delete_object!(vis,scope,"ellipse_contour")
	end

	# Animate the scene
	anim = MeshCat.Animation(anim.clips, opts.framerate)
	scene_animation(solver, scenario, vis, anim, opts)
    MeshCat.setanimation!(vis, anim)
	return vis, anim
end

function build_actors(vis::Visualizer, scenario::Scenario{T}, scope::Symbol,
	opts::AnimationOptions{T}=AnimationOptions{T}()) where T
	scale = opts.scale
	α_fading = opts.α_fading

	actors_radii = scenario.actors_radii*scale
	actors_types = scenario.actors_types

	car_offset = compose(compose(compose(
			Translation(-0.03, 0.0, 0.31),
			LinearMap(AngleAxis(-0.012+pi/200, 0, 1, 0))),
			LinearMap(AngleAxis(pi/2, 0, 0, 1))),
			LinearMap(AngleAxis(pi/2, 1, 0, 0)))
	ped_offset = LinearMap(AngleAxis(pi/2, 0, 0, 1))
	cylinder_height_car = 0.01*scale
	cylinder_height_ped = 0.035*scale

	pkg_path = joinpath(dirname(@__FILE__), "../../")

	path = String(scope)*"/actors/actor_bundle_"
	for i = 1:length(actors_radii)
		col_name = opts.color_set[i]
		col_val = animation_colors(i,opts)
		if actors_types[i] == :car
			cylinder_height = cylinder_height_car
			obj_path = joinpath(pkg_path, "resources/object/car_geometry.obj")
			mtl_path = joinpath(pkg_path, "resources/material/"*col_name*".mtl")
			actor = ModifiedMeshFileObject(obj_path, mtl_path, scale=0.106)
			setobject!(vis[path*"$i/actor"], actor)
			settransform!(vis[path*"$i/actor"], car_offset)
		elseif actors_types[i] == :pedestrian
			cylinder_height = cylinder_height_ped
			# actor = load(joinpath(pkg_path, "resources/object/pedestrian_geometry.obj"), GLUVMesh)
			# actor.vertices .*= 0.001/15 * scale
			# actor_material = MeshPhongMaterial(color=RGBA(0., 0., 0., 1.0))
			obj_path = joinpath(pkg_path, "resources/object/pedestrian_geometry.obj")
			mtl_path = joinpath(pkg_path, "resources/material/gray.mtl")
			actor = ModifiedMeshFileObject(obj_path, mtl_path, scale=0.001/15 * scale)
			setobject!(vis[path*"$i/actor"], actor)
			# setobject!(vis[path*"$i/actor"], actor, actor_material)
			settransform!(vis[path*"$i/actor"], ped_offset)
		end
		# Add collision avoidance cylinders
		collision_cylinder = Cylinder(Point(0.0,0.0,0.0),
		Point(0.0,0.0,cylinder_height), actors_radii[i])
		cylinder_material = MeshPhongMaterial(color=RGBA(col_val..., α_fading))
		setobject!(vis[path*"$i/col_cyl"], collision_cylinder, cylinder_material)
	end
	return nothing
end

function build_trajectory(vis::Visualizer, solver::TO.AbstractSolver{T},
	opts::AnimationOptions{T}=AnimationOptions{T}()) where {T}
	build_trajectory(vis,solver,TO.states(solver),opts)
end

function build_trajectory(vis::Visualizer, solver::TO.AbstractSolver{T},
	X::Vector{K}, opts::AnimationOptions{T}=AnimationOptions{T}()) where {T,K}

	dt = opts.dt
	scale = opts.scale
	α_fading = opts.α_fading

	cyl_height = 0.005*scale
	cyl_radius = 0.015*scale

	model = TO.get_model(solver)
	p = model.p
	px = model.px
	path = String(solver_scope(solver, opts))*"/trajectory"
	for k = 1:length(X)
		for i = 1:p
			col_val = animation_colors(i,opts)
			cyl_pos = Translation(scale*X[k][px[i]]..., 0)
			aug_cyl_height = cyl_height + i*0.00 + k*0.01/length(X)*scale
			cyl = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,aug_cyl_height), cyl_radius)
			cyl_material = MeshPhongMaterial(color=RGBA(col_val..., α_fading))
			setobject!(vis[path*"/cyl_$i/step"*"_$k"], cyl, cyl_material)
			settransform!(vis[path*"/cyl_$i/step"*"_$k"], cyl_pos)
		end
	end
	return nothing
end

function build_open_loop_plan(vis::Visualizer, solver::TO.AbstractSolver{T},
	opts::AnimationOptions{T}=AnimationOptions{T}()) where {T,K}

	horizon = opts.horizon
	dt = opts.dt
	scale = opts.scale
	α_fading = opts.α_fading

	cyl_height = 0.007*scale
	cyl_radius = 0.010*scale

	model = TO.get_model(solver)
	p = model.p
	path = String(solver_scope(solver, opts))*"/open_loop_plan"

	for j = 1:horizon
		for i = 1:p
			col_val = animation_colors(i,opts)
			aug_cyl_height = cyl_height # + i*0.00 + k*0.01/M*scale
			cyl = Cylinder(Point(0.0,0.0,0.0),
				Point(0.0,0.0,aug_cyl_height), cyl_radius)
			cyl_material = MeshPhongMaterial(color=RGBA(col_val..., α_fading))
			setobject!(vis[path*"/cyl_$i/stage_$j"], cyl, cyl_material)
		end
	end
	return nothing
end

function build_sampled_plan(vis::Visualizer, solver::TO.AbstractSolver{T},
	opts::AnimationOptions{T}=AnimationOptions{T}()) where {T,K}

	horizon = opts.horizon
	dt = opts.dt
	scale = opts.scale
	α_fading = opts.α_fading
	col_val = animation_colors("red")

	cyl_height = 0.020*scale
	cyl_radius = 0.005*scale

	model = TO.get_model(solver)
	p = model.p
	num_samples = 2*length(solver.θ_inds)+1
	path = String(solver_scope(solver, opts))*"/sampled_plan"
	for l = 1:num_samples
		for j = 1:horizon
			for i = 1:p
				aug_cyl_height = cyl_height# + i*0.00 + k*0.01/M*scale
				cyl = Cylinder(Point(0.0,0.0,0.0),
					Point(0.0,0.0,aug_cyl_height), cyl_radius)
				cyl_material = MeshPhongMaterial(color=RGBA(col_val..., α_fading))
				setobject!(vis[path*"/cyl_$i/stage_$j/sample_$l"], cyl, cyl_material)
			end
		end
	end
	return nothing
end

# function build_ellipse_contour(vis::Visualizer, solver::TO.AbstractSolver{T},
# 	opts::AnimationOptions{T}=AnimationOptions{T}()) where {T,K}
#
# 	scale = opts.scale
# 	α_fading = opts.α_fading
# 	col_val = animation_colors("yellow")
#
# 	cyl_height = 0.050*scale
# 	cyl_radius = 0.005*scale
#
# 	model = TO.get_model(solver)
# 	p = model.p
# 	path = String(solver_scope(solver, opts))*"/ellipse_contour"
# 	for l = 1:opts.num_contour_points
# 		for i in setdiff(1:p, solver.i_av)
# 			cyl = Cylinder(Point(0.0,0.0,0.0),
# 				Point(0.0,0.0,cyl_height), cyl_radius)
# 			cyl_material = MeshPhongMaterial(color=RGBA(col_val..., α_fading))
# 			setobject!(vis[path*"/cyl_$i/sample_$l"], cyl, cyl_material)
# 		end
# 	end
# 	return nothing
# end

function build_ellipse_contour(vis::Visualizer, solver::TO.AbstractSolver{T},
	scenario::Scenario{T}, opts::AnimationOptions{T}=AnimationOptions{T}()) where {T,K}

	scale = opts.scale
	α_fading = opts.α_fading
	horizon = opts.horizon
	algo_type = solver.opts.algorithm_type
	start_ind, end_ind, α = time_resampling([z.dt for z in solver.stats.traj], dt=opts.dt)
	model = TO.get_model(solver)
	p = model.p
	px = model.px
	path = String(solver_scope(solver, opts))*"/ellipse_contour"

	col_val = animation_colors("yellow")
	ellipse_thickness = 0.001*scale

	for q = 1:length(α)
		for k = 1:horizon
			ellipse_height = 0.001*scale*(1+k/5)
			col_val_k = col_val + [-k*0.04, 0, 0]
			obj_material = MeshPhongMaterial(color=RGBA(col_val_k..., 0.5))#α_fading/10))
			start_ell = solver.stats.traj_ellip[start_ind[q]][k]
			end_ell = solver.stats.traj_ellip[end_ind[q]][k]
			c = (1-α[q]) * start_ell.mean + α[q] * end_ell.mean
			# Σ = (1-α[q]) * start_ell.cov + α[q] * end_ell.cov###############################################
			Σ = (1-1.) * start_ell.cov + 1. * end_ell.cov
			ell = Ellipsoid(c, Σ)
			for i in setdiff(1:p, solver.i_av)
				inds = px[i]
				radius = scenario.actors_radii[i]
				ell_i = Ellipsoid(SVector{length(px[i])}(ell.mean[inds]), ell.cov[inds, inds])
				obj = build_ellipse_object(ell_i, radius; scale=opts.scale, height=ellipse_height, thickness=ellipse_thickness)
				setobject!(vis[path*"/ell_$i/step_$k/time_$q"], obj, obj_material)
			end
		end
	end
	return nothing
end

function scene_animation(solver::TO.AbstractSolver, scenario::Scenario{T}, vis::Visualizer,
	anim::MeshCat.Animation, opts::AnimationOptions{T}=AnimationOptions{T}())  where T

	actor_tracking = opts.actor_tracking
	tracking_ind = opts.tracking_ind
	horizon = opts.horizon
	dt = opts.dt
	scale = opts.scale
	camera_offset = opts.camera_offset
	camera_mvt = opts.camera_mvt
	scope = solver_scope(solver, opts)

	n,m,N = size(solver)
	model = TO.get_model(solver)
	p = model.p
	px = model.px
	offset = camera_offset ? -0.80*scale : 0.00

	# Animate the scene
	birdseye_trans = Translation(0.0, 0.0, 1.9*scale)
    birdseye_rot = compose(
		LinearMap(AngleAxis(-pi/2, 0, 0, 1)),
		LinearMap(AngleAxis(-0.397*pi, 0, 1, 0)))


    # Compute actor transformations
    actor_translations, actor_rotations = actor_transformations(solver, opts)
	if occursin("UKFGamesSolver", string(scope))
		start_ind, end_ind, α = time_resampling([z.dt for z in solver.stats.traj], dt=dt)
	end

	# We get rid of the last frame if the final state constraints are relaxed
	actor_path = String(scope)*"/actors/actor_bundle"
	trajectory_path = String(scope)*"/trajectory"
	open_loop_plan_path = String(scope)*"/open_loop_plan"
	sampled_plan_path = String(scope)*"/sampled_plan"
	ellipse_contour_path = String(scope)*"/ellipse_contour"
	roadway_path = "roadway"

	for k = 1:length(actor_translations)-1
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

			if occursin("UKFGamesSolver", string(scope))
				num_samples = 2*length(solver.θ_inds)+1
				for j=1:horizon
					x_start = TO.state(solver.stats.traj_pred[start_ind[k]][j])
					x_end = TO.state(solver.stats.traj_pred[end_ind[k]][j])
					x = x_start + α[k]*(x_end-x_start)
					for i=1:p
						cyl_pos = compose(roadway_translation, Translation(scale*x[px[i]]..., 0))
						settransform!(vis[open_loop_plan_path*"/cyl_$i/stage_$j"], cyl_pos)
					end
				end
				if opts.display_sampled_plan
					if k == 1
						for j=1:horizon
							x_start = TO.state(solver.stats.traj_pred[start_ind[k]][j])
							x_end = TO.state(solver.stats.traj_pred[end_ind[k]][j])
							x = x_start + α[k]*(x_end-x_start)
							for i=1:p
								cyl_pos = compose(roadway_translation, Translation(scale*x[px[i]]..., 0))
								for l = 1:num_samples
									settransform!(vis[sampled_plan_path*"/cyl_$i/stage_$j/sample_$l"], cyl_pos)
								end
							end
						end
					else
						for l = 1:num_samples
							for j=1:horizon
								x_start = TO.state(solver.stats.traj_sampled[start_ind[k]+1][l][j])
								x_end = TO.state(solver.stats.traj_sampled[end_ind[k]+1][l][j])
								x = x_start + α[k]*(x_end-x_start)
								for i=1:p
									cyl_pos = compose(roadway_translation, Translation(scale*x[px[i]]..., 0))
									settransform!(vis[sampled_plan_path*"/cyl_$i/stage_$j/sample_$l"], cyl_pos)
								end
							end
						end
					end
				end
				# if opts.display_ellipse_contour &&
				# 	solver.opts.algorithm_type == :algames &&
				# 	solver.opts.ellipsoid_collision
				#
				# 	if k > 2
				# 		for i in setdiff(1:p, solver.i_av)
				# 			cont_start = ellipse_trajectory_contour(
				# 				solver.stats.traj_ellip[start_ind[k]],
				# 				px[i],
				# 				scenario.actors_radii[i];
				# 				P=opts.num_contour_points,
				# 				M=50)
				# 			cont_end = ellipse_trajectory_contour(
				# 				solver.stats.traj_ellip[end_ind[k]],
				# 				px[i],
				# 				scenario.actors_radii[i];
				# 				P=opts.num_contour_points,
				# 				M=50)
				# 			for l = 1:opts.num_contour_points
				# 				cont = cont_start[l,:] + α[k]*(cont_end[l,:]-cont_start[l,:])
				# 				@show cont
				# 				cyl_pos = compose(roadway_translation, Translation(scale*cont..., 0))
				# 				settransform!(vis[ellipse_contour_path*"/cyl_$i/sample_$l"], cyl_pos)
				# 			end
				# 		end
				# 	end
				# end
				if opts.display_ellipse_contour &&
					solver.opts.algorithm_type in [:ground_truth, :line_prediction, :algames, :frozen_learning] &&
					solver.opts.ellipsoid_collision
					for k = 1:length(α)
						for j = 1:horizon
							for i in setdiff(1:p, solver.i_av)
								setvisible!(vis[ellipse_contour_path*"/ell_$i/step_$j/time_$k"], false)
							end
						end
					end

					if k > 2
						for j = 1:horizon
							start_ell = solver.stats.traj_ellip[start_ind[k]][j]
							end_ell = solver.stats.traj_ellip[end_ind[k]][j]
							c = (1-α[k]) * start_ell.mean + α[k] * end_ell.mean
							Σ = (1-α[k]) * start_ell.cov + α[k] * end_ell.cov
							ell = Ellipsoid(c, Σ)
							for i in setdiff(1:p, solver.i_av)
								inds = px[i]
								radius = scenario.actors_radii[i]
								ell_i = Ellipsoid(SVector{length(px[i])}(ell.mean[inds]), ell.cov[inds, inds])
								transform = ellipse_transformation(ell_i, radius, scale=opts.scale)
								ell_pos = compose(roadway_translation, transform)
								settransform!(vis[ellipse_contour_path*"/ell_$i/step_$j"], ell_pos)
								setvisible!(vis[ellipse_contour_path*"/ell_$i/step_$j/time_$k"], true)
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
			camera_transformation = compose(
				birdseye_trans,
				birdseye_rot)
            ##
			# camera_transformation = compose(compose(
			# 	Translation(-4.0, 0.0, -k*0.6),
			# 	compose(
			# 	birdseye_trans,
			# 	birdseye_rot)),Translation(0.0, -k*0.3, -5+k*0.6))
			# ##
			settransform!(vis["/Cameras/default"], camera_transformation)

			# ##
			# setprop!(
			# 	# vis["/Cameras/default"],
			# 	vis["/Cameras/default/rotated/<object>"],
			# 	"position",
			# 	[-15.0 + 1.0*k, 0., 0.])
			# ##
			##
			# setprop!(
			# 	# vis["/Cameras/default"],
			# 	vis["/Cameras/default/rotated/<object>"],
			# 	"position",
			# 	[-15., -5., +10.])
			# ##
        end
	end
	setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 1.0)
	setprop!(vis["/Lights/DirectionalLight/<object>"], "intensity", 1.2)
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 3.0)
	return nothing
end

function actor_transformations(solver::TO.AbstractSolver{T},
	opts::AnimationOptions{T}=AnimationOptions{T}()) where {T}
	actor_transformations(solver,TO.states(solver),opts)
end

function actor_transformations(solver::TO.AbstractSolver{T}, X::Vector{K},
	opts::AnimationOptions{T}=AnimationOptions{T}()) where {T,K}

	scale = opts.scale
	n,m,N = size(solver)
	px = TO.get_model(solver).px
	p = length(px)
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

function time_resampling(dts::Vector{T}; dt::T=5e-1) where {T}
	Δt = sum(dts)
	sum_dt = cumsum([0.;dts])
	M = Int(floor(Δt/dt))
	# M = Int(floor(Δt/dt)-1)
	t = 0.
	start_ind = zeros(Int, M)
	end_ind = zeros(Int, M)
	α = zeros(T, M)
	for k = 1:M
		start_ind[k] = findlast(x -> x <= t, sum_dt)
		end_ind[k] = min(start_ind[k]+1, length(sum_dt)-1)
		α[k] = (t - sum_dt[start_ind[k]])/(sum_dt[end_ind[k]] - sum_dt[start_ind[k]])
		t += dt
	end
	return start_ind, end_ind, α
end

function state_resampling(traj::Vector{K}; dt::T=5e-1) where {T,K}
	n = length(traj[1]._x)
	start_ind, end_ind, α = time_resampling([z.dt for z in traj], dt=dt)
	M = length(α)
	X = Vector{SVector{n,T}}(undef, M)
	for k = 1:M
		x_start = TO.state(traj[start_ind[k]])
		x_end = TO.state(traj[end_ind[k]])
		X[k] = x_start + α[k]*(x_end-x_start)
	end
	return X
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


function build_ellipse_object(ell::Ellipsoid{n,T}, radius::T; scale::T=0.1, s::T=5.991, height::T=0.0, thickness::T=0.05) where {n,T}
	@assert n == 2
	c = ell.mean
	Σ = ell.cov
	inflation = radius*ones(length(c))
	Splot = inv(sqrt(s*Σ) + Diagonal(inflation))
	λ1, λ2 = eigvals(Splot)

	obj = HyperEllipsoid{3,T}(
		Point(0., 0., height),
		Vec(scale*1/λ1, scale*1/λ2, thickness)
		)
	return obj
end

function ellipse_transformation(ell::Ellipsoid{n,T}, radius::T; scale::T=0.1, s::T=5.991) where {n,T}
	@assert n == 2
	c = ell.mean
	Σ = ell.cov
	inflation = radius*ones(length(c))
	Splot = inv(sqrt(s*Σ) + Diagonal(inflation))
	e1 = eigvecs(Splot)[:,1]
	θ1 = atan(e1[2], e1[1])

	rot = LinearMap(AngleAxis(θ1, 0, 0, 1))
	trans = Translation(scale .* c..., 0.)
	transform = compose(trans, rot)
	return transform
end
