################################################################################
# Utils
################################################################################

# rotation matrix rotating unit vector a onto unit vector b
function rot(a, b)
	v = cross(a, b)
	s = sqrt(transpose(v) * v)
	c = transpose(a) * b

	R = Diagonal(@SVector ones(3)) + skew(v) + 1.0 / (1.0 + c) * skew(v) * skew(v)
end

function skew(x)
	SMatrix{3,3}([0.0 -x[3] x[2];
	               x[3] 0.0 -x[1];
				   -x[2] x[1] 0.0])
end

function default_background!(vis)
    setvisible!(vis["/Background"], true)
    setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
    setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
    setvisible!(vis["/Axes"], false)
end

function get_material(;α=1.0)
    orange_col = [1,153/255,51/255]
    blue_col = [51/255,1,1]
    black_col = [0,0,0]
    orange_mat = MeshPhongMaterial(color=MeshCat.RGBA(orange_col...,α))
    blue_mat = MeshPhongMaterial(color=MeshCat.RGBA(blue_col...,α))
    black_mat = MeshPhongMaterial(color=MeshCat.RGBA(black_col...,α))
    return orange_mat, blue_mat, black_mat
end

################################################################################
# Static Object
################################################################################

function build_traj!(vis::Visualizer, model::AbstractGameModel, traj::Traj; name::Symbol=:Traj, α=0.8, r=0.02)
	default_background!(vis)
	r = convert(Float32, r)
    p = model.p
    pz = model.pz
	d = dim(model)
	orange_mat, blue_mat, black_mat = get_material(;α=α)

    for t in eachindex(traj)
        s = TrajectoryOptimization.state(traj[t])
        for i = 1:p
            x = [s[pz[i][1:d]]; zeros(3-d)]
            setobject!(vis[name]["player$i"]["t$t"], Sphere(Point3f0(0.0), r), blue_mat)
            settransform!(vis[name]["player$i"]["t$t"], MeshCat.Translation(x...))
        end
    end
    return nothing
end

function build_xf!(vis::Visualizer, model::AbstractGameModel, xf::AbstractVector; name::Symbol=:Xf, α=0.8)
	default_background!(vis)
    p = model.p
    pz = model.pz
	d = dim(model)
	orange_mat, blue_mat, black_mat = get_material(;α=α)

    for i = 1:p
        x = [xf[i][1:d]; zeros(3-d)]
        setobject!(vis[name]["player$i"], Sphere(Point3f0(0.0), 0.05), orange_mat)
        settransform!(vis[name]["player$i"], MeshCat.Translation(x...))
    end
    return nothing
end

################################################################################
# Animations
################################################################################

function animate_robot!(vis::Visualizer, anim::MeshCat.Animation, model::AbstractGameModel,
		s::AbstractVector; name::Symbol=:Robot)
	default_background!(vis)
	for t in 1:length(s)
		MeshCat.atframe(anim, t) do
			set_robot!(vis, model, s[t], name=name)
		end
	end
	setanimation!(vis, anim)
	return nothing
end

function visualize_robot!(vis::Visualizer, model::AbstractGameModel, s::AbstractVector;
		h=0.01, α=1.0,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=:Robot)

	build_robot!(vis, model, name=name, α=α)
	animate_robot!(vis, anim, model, s, name=name)
	return anim
end

function visualize_robot!(vis::Visualizer, model::AbstractGameModel, traj::Traj;
		sample=max(1, Int(floor(traj[1].dt*length(traj) / 100))), h=traj[1].dt*sample, α=1.0,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=:Robot)

	visualize_robot!(vis, model, TrajectoryOptimization.state.(traj); anim=anim, name=name, h=h, α=α)
	return anim
end
