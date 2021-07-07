# function build_wall!(vis::Visualizer, walls::AbstractVector{<:AbstractWall}; name::Symbol=:Env, α=0.2)
#     for i in eachindex(walls)
#         build_wall!(vis, walls[i], name=name, id=i, α=α)
#     end
#     return nothing
# end
#
# function build_cylinder!(vis::Visualizer, walls::AbstractVector{<:AbstractWall}; radius=0.0, name::Symbol=:Env, α=0.2)
#     for i in eachindex(walls)
#         build_cylinder!(vis, walls[i], radius=radius, name=name, id=i, α=α)
#     end
#     return nothing
# end
#
# function build_wall!(vis::Visualizer, wall::Wall3D; name::Symbol=:Env, id::Int, α=0.2)
#     p1 = wall.p1
#     p2 = wall.p2
#     p3 = wall.p3
#     p4 = p3 + (p1 - p2)
#     l1 = norm(p1 - p2)
#     l2 = norm(p3 - p2)
#     v1 = (p1 - p2)/l1
#     v2 = (p3 - p2)/l2
#     orange_mat, blue_mat, black_mat = get_material(;α=α)
#
#
#     setobject!(vis[name]["wall$id"][:slab],
#         MeshCat.HyperRectangle(MeshCat.Vec(0.0, 0.0, 0.0), MeshCat.Vec(l1, l2, 0.0)),
#         black_mat)
#     for j = 1:3
#         setobject!(vis[name]["wall$id"]["p$j"],
#             GeometryBasics.Sphere(MeshCat.Point3f0(0.0), 0.05),
#             black_mat)
#     end
#
#     rt = MeshCat.LinearMap([v1 v2 cross(v1, v2)])
#     settransform!(vis[name]["wall$id"][:slab], MeshCat.compose(MeshCat.Translation(p2...), rt))
#     settransform!(vis[name]["wall$id"]["p1"], MeshCat.Translation(p1...))
#     settransform!(vis[name]["wall$id"]["p2"], MeshCat.Translation(p2...))
#     settransform!(vis[name]["wall$id"]["p3"], MeshCat.Translation(p3...))
#     return nothing
# end
#
# function build_cylinder!(vis::Visualizer, wall::CylinderWall; radius=0.0, name::Symbol=:Env, id::Int, α=0.2)
#     p = wall.p
#     v = wall.v
#     v_ = Int.([v == :x, v == :y, v == :z])
#     l = convert(Float32, wall.l)
#     r = convert(Float32, wall.r-radius)
#     orange_mat, blue_mat, black_mat = get_material(;α=α)
#
#     setobject!(vis[name]["cylinder$id"][:body],
#         MeshCat.Cylinder(MeshCat.Point3f0(0.0), MeshCat.Point3f0((l*v_)...), r),
#         black_mat)
#
#     # MeshCat.HyperRectangle(MeshCat.Vec(0.0, 0.0, 0.0),
#     #     MeshCat.Vec(l1, l2, 0.0)), MeshPhongMaterial(color = MeshCat.RGBA(0.5, 0.5, 0.5, α)))
#     setobject!(vis[name]["cylinder$id"]["p"], GeometryBasics.Sphere(MeshCat.Point3f0(0.0), 0.05), black_mat)
#
#     # rt = MeshCat.LinearMap([v1 v2 cross(v1, v2)])
#     # settransform!(vis[name]["cylinder$id"][:body], MeshCat.compose(MeshCat.Translation(p...), rt))
#     settransform!(vis[name]["cylinder$id"][:body], MeshCat.Translation(p...))
#     settransform!(vis[name]["cylinder$id"]["p"], MeshCat.Translation(p...))
#     return nothing
# end
