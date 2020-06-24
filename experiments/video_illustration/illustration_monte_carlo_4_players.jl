using Statistics
using Random
using Plots
using Logging
using ALGAMES
using PGFPlotsX
const AG = ALGAMES

include("algames_solver_4_players.jl")
include("approx_solver_4_players.jl")
include("const_vel_solver_4_players.jl")
include("mpc_vect_solver.jl")
include("mpc_vect_methods.jl")

# Algames Solver
ramp_merging_4_players_unicycle_penalty_scenario
algames_solver
approx_solvers
const_vel_solvers

dxf = @SVector [
               0.60, 0.60, 0.70, 0.70,
			   0.00, 0.00, 0.00, 0.00,
			   0.00, 0.00, 0.00, 0.00,
			   0.00, 0.00, 0.00, 0.00,
			   ]
state_noise = 0. * SVector{n}([ #######################################################################
	0.008,   0.008,   0.008,   0.008,
	0.008,   0.008,   0.008,   0.008,
	2*pi/72, 2*pi/72, 2*pi/72, 2*pi/72,
	0.03,    0.03,    0.03,    0.03]) #+-50cm, +-50cm, +-25deg, +-12.5% per second
timeout = 0.400
opts_mpc = MPCGamesSolverOptions{n,T}(
	# live_plotting=:on,
	iterations=2000,
	N_mpc=100,
	mpc_tf=6.0,
	min_δt=0.005,
	max_δt=timeout,
    dxf=dxf,
	noise=state_noise)

algames_solver.opts.timeout = timeout
algames_mpc_solver = MPCGamesSolver(algames_solver, opts_mpc)
reset!(algames_mpc_solver, reset_type=:full)
solve!(algames_mpc_solver; wait=false)
resample!(algames_mpc_solver)
algames_mpc_solver.stats.failure_status

using MeshCat
algames_vis = MeshCat.Visualizer()
algames_anim = MeshCat.Animation()
open(algames_vis)
sleep(1.0)
algames_vis, algames_anim = AG.animation(algames_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
	vis=algames_vis, anim=algames_anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true,
	camera_offset=false)


# approx_mpc_solver = MPCVectGamesSolver11(approx_solvers, opts_mpc)
# reset!(approx_mpc_solver, reset_type=:full)
# solve!(approx_mpc_solver; wait=false)
# resample!(approx_mpc_solver)
# approx_mpc_solver.stats.failure_status

const_vel_mpc_solver = MPCVectGamesSolver11(const_vel_solvers, opts_mpc)
reset!(const_vel_mpc_solver, reset_type=:full)
solve!(const_vel_mpc_solver; wait=false)
resample!(const_vel_mpc_solver)
const_vel_mpc_solver.stats.failure_status


using MeshCat
const_vel_vis = MeshCat.Visualizer()
const_vel_anim = MeshCat.Animation()
open(const_vel_vis)
sleep(1.0)
const_vel_vis, const_vel_anim = animation(const_vel_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
	vis=const_vel_vis, anim=const_vel_anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true,
	display_open_loop_plan=true,
	camera_offset=false)


maximum(algames_mpc_solver.stats.cmax_collision)
maximum(algames_mpc_solver.stats.cmax_boundary)
maximum(algames_mpc_solver.stats.cmax_bound)

maximum(approx_mpc_solver.stats.cmax_collision)
maximum(approx_mpc_solver.stats.cmax_boundary)
maximum(approx_mpc_solver.stats.cmax_bound)

maximum(const_vel_mpc_solver.stats.cmax_collision)
maximum(const_vel_mpc_solver.stats.cmax_boundary)
maximum(const_vel_mpc_solver.stats.cmax_bound)



mpc_state_noise = 1.0*SVector{n}([
	0.20,    0.20,    0.20,    0.30,  # 5  m
	0.02,    0.02,    0.02,    0.02,  # 50 cm
	pi/30,   pi/30,   pi/30,   pi/30, # 6  deg
	0.01,    0.01,    0.01,    0.01]) # 1  m/s
iter = 100
algames_failure, algames_rank = monte_carlo_analysis(
	algames_solver,
	opts_mpc,
	mpc_state_noise,
	iter)

# approx_failure, approx_rank = monte_carlo_analysis(
# 	approx_solvers,
# 	opts_mpc,
# 	mpc_state_noise,
# 	iter)

const_vel_failure, const_vel_rank = monte_carlo_analysis(
	const_vel_solvers,
	opts_mpc,
	mpc_state_noise,
	iter)

const_vel_failure_save = deepcopy(const_vel_failure)
const_vel_rank_save = deepcopy(const_vel_rank)

visualize_rank(algames_rank[.!algames_failure])
visualize_rank(const_vel_rank[.!const_vel_failure])


algames_vis=AG.Visualizer()
algames_anim=AG.MeshCat.Animation()
open(algames_vis)
# Execute this line after the MeshCat tab is open
algames_vis, algames_anim = AG.animation(algames_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
	vis=algames_vis, anim=algames_anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true,
	camera_offset=false)


# approx_vis=AG.Visualizer()
# approx_anim=AG.MeshCat.Animation()
# open(approx_vis)
# Execute this line after the MeshCat tab is open
approx_vis, approx_anim = animation(approx_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
	vis=approx_vis, anim=approx_anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true,
	camera_offset=false)

const_vel_vis=AG.Visualizer()
const_vel_anim=AG.MeshCat.Animation()
open(const_vel_vis)
const_vel_vis, const_vel_anim = animation(const_vel_mpc_solver, ramp_merging_4_players_unicycle_penalty_scenario;
	vis=const_vel_vis, anim=const_vel_anim,
	open_vis=false,
	display_actors=true,
	display_trajectory=true,
	camera_offset=false)



iter = algames_mpc_solver.stats.iterations
mean_solve_time = mean(algames_mpc_solver.stats.solve_time[1:iter])
update_freq = algames_mpc_solver.stats.iterations/algames_mpc_solver.stats.time
std_solve_time = std(algames_mpc_solver.stats.solve_time)
largest_solve_time = maximum(algames_mpc_solver.stats.solve_time)

# iter = const_vel_mpc_solver.stats.iterations
# mean_solve_time = mean(const_vel_mpc_solver.stats.solve_time[1:iter])
# update_freq = const_vel_mpc_solver.stats.iterations/mpc_solver.stats.time
# std_solve_time = std(const_vel_mpc_solver.stats.solve_time)
# largest_solve_time = maximum(const_vel_mpc_solver.stats.solve_time)


# algames_mpc_solver.opts.N_mpc = 500
# approx_mpc_solver.opts.N_mpc = 500
# const_vel_mpc_solver.opts.N_mpc = 500
#
visualize_control(algames_mpc_solver)
visualize_control(const_vel_mpc_solver)
visualize_state(algames_mpc_solver)
visualize_state(const_vel_mpc_solver)
visualize_control(const_vel_solvers[1])
visualize_control(const_vel_solvers[2])
visualize_control(const_vel_solvers[3])
visualize_control(const_vel_solvers[4])


minimum([state[end] for state in TO.states(const_vel_mpc_solver)[1:end-2]])

# using FileIO
# save("const_vel_mpc_solver.Z.jld2", "const_vel_mpc_solver.Z", const_vel_mpc_solver.Z)
# save("algames_mpc_solver.Z.jld2", "algames_mpc_solver.Z", algames_mpc_solver.Z)
#
# const_vel_mpc_solver.Z

using PGFPlotsX
function visualize_latex_speed(algames_mpc_solver::MPCGamesSolver{T},
	const_vel_mpc_solver::MPCVectGamesSolver11{T}; save=false) where T

	N_final = Int(floor(0.75*algames_mpc_solver.opts.N_mpc))
	algames_states = TO.states(algames_mpc_solver)[1:N_final]
	const_vel_states = TO.states(const_vel_mpc_solver)[1:N_final]

	algames_speeds = [s[end] for s in algames_states]
	const_vel_speeds = [s[end] for s in const_vel_states]

    max_vel = max(maximum(algames_speeds), maximum(const_vel_speeds))
    algames_speeds ./= max_vel
    const_vel_speeds ./= max_vel


	time = range(0, length=length(algames_speeds), stop=1.0)
    x = range(-1; stop = 1, length = 51) # so that it contains 1/0
    axis = @pgf Axis(
        {
            "legend pos=south east",
            ymajorgrids,
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
            "x label style={at={(axis description cs:0.5,-0.20)},anchor=north}",
            "y label style={at={(axis description cs:-0.10,0.5)},rotate=0,anchor=south}",
            "xlabel={Scaled Time}",
            "ylabel={Scaled Velocity}",
            xmajorgrids = false,
            xmin = 0.00,   xmax = 1.00,
            ymin = 0.00, #  ymax = 1.00,

        },
        Plot(
            {
            "thick",
            "orange",
            no_marks,
            },
            Coordinates(time, algames_speeds)
        ),
        LegendEntry("MPC ALGAMES"),
        Plot(
            {
            "thick",
            "blue",
            no_marks,
            },
            Coordinates(time, const_vel_speeds)
        ),
        LegendEntry("MPC Baseline"),
    )
    pgfsave("plots/tikz/velocity.tikz",
        axis; include_preamble=false, dpi = 600)
    return axis
end

function visualize_rank(vec_rank)
	plt = plot()
    histogram!(vec_rank,
        bins=10,
        color=:green,
        linewidth=1.0,
        legend=:bottomright,
        title="Vehicle Ordering",
        xlabel="Merging Vehicle's Place",
        ylabel="Occurences",
        label="")
    display(plt)
	return nothing
end


visualize_latex_speed(algames_mpc_solver, const_vel_mpc_solver, save=true)

function visualize_latex_rank(algames_rank, const_vel_rank; save=false) where T
	labels = ["2nd Place", "3rd Place", "4th Place"]
	p = 4
	algames_occ = zeros(Int, p)
	const_vel_occ = zeros(Int, p)
	for k = 1:p
		algames_occ[k] = length(findall(x->x==k, algames_rank))
		const_vel_occ[k] = length(findall(x->x==k, const_vel_rank))
	end
	algames_coords = [(algames_occ[i+1], labels[i]) for i=length(labels):-1:1]
	# algames_coords = [(labels[i], algames_occ[i+1]) for i=1:length(labels)]
	const_vel_coords = [(const_vel_occ[i+1], labels[i]) for i=length(labels):-1:1]
	# const_vel_coords = [(labels[i], const_vel_occ[i+1]) for i=1:length(labels)]
	axis = @pgf Axis(
	    {
	        xbar,
	        enlargelimits = 0.15,
	        legend_style =
	        {
	            at = Coordinate(0.5, -0.15),
	            anchor = "north",
	            legend_columns = -1
	        },
	        # ylabel = raw"\#participants",
	        symbolic_y_coords=labels,
			ytick = "data",
			xmin = 0.00,
	        nodes_near_coords,
	        nodes_near_coords_align={horizontal},
	    },
		Plot({"draw=black", "fill=orange!50"},
			Coordinates(const_vel_coords)),
		Plot({"draw=black", "fill=blue!25"},
			Coordinates(algames_coords)),
	    Legend(["Baseline", "ALGAMES"])
	)
	pgfsave("plots/tikz/mpc_rank2.tikz",
		axis; include_preamble=false, dpi = 600)
	return axis
end

visualize_latex_rank(algames_rank, const_vel_rank)
using StatsBase
