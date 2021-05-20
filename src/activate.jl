# Run tests locally
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.test("Algames")


Pkg.activate(joinpath(@__DIR__, "../test"))

# Pkg.add(Pkg.PackageSpec(;name="RobotDynamics", version="0.3.1"))
# Pkg.add(Pkg.PackageSpec(;name="TrajectoryOptimization", version="0.4.1"))
# Pkg.add(Pkg.PackageSpec(;name="Altro", version="0.3.0"))

# Pkg.add("Plots")
# using Plots
# Pkg.remove("Plots")

# Pkg.add("MeshCat")
# using MeshCat
# Pkg.remove("MeshCat")

# Pkg.add("GeometryBasics")
# using GeometryBasics
# Pkg.remove("GeometryBasics")
