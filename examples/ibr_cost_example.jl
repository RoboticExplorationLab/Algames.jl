################################################################################
# Iterative Best Response Example
################################################################################
# using Algames
using Plots
using StaticArrays
using LinearAlgebra

function ibr_experiment(p::Int)
	T = Float64

	# Define the dynamics of the system
	d = 2
	model = DoubleIntegratorGame(p=p, d=d) # game with 3 players with unicycle dynamics
	n = model.n
	m = model.m

	# Define the horizon of the problem
	N = 20 # N time steps
	dt = 0.1 # each step lasts 0.1 second
	probsize = ProblemSize(N,model) # Structure holding the relevant sizes of the problem

	# Define the objective of each player
	# We use a LQR cost
	Q = [Diagonal(10*ones(SVector{model.ni[i],T})) for i=1:p] # Quadratic state cost
	R = [Diagonal(0.01*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
	# Desrired state
	xf = [SVector{model.ni[1],T}([ 2, 0.2*i,0,0]) for i = 0:p-1]
	# Desired control
	uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
	# Objectives of the game
	game_obj = GameObjective(Q,R,xf,uf,N,model)
	radius = 1.0*ones(p)
	μ = 15.0*ones(p)
	add_collision_cost!(game_obj, radius, μ)

	# Define the constraints that each player must respect
	game_con = GameConstraintValues(probsize)

	# Define the initial state of the system
	x0 = [[ 0.0,  0.2*i,  0.0,  0.0] for i=0:p-1]
	x0 = SVector{model.n,T}(reshape(vcat([x0[i]' for i=1:p]...), model.n))

	# Define the Options of the solver
	opts_alg = Options(Δ_min=1e-9, inner_iter=2, inner_print=false, ϵ_opt=1e-9)
	opts_ibr = Options(dual_reset=false, outer_iter=1, inner_iter=1, Δ_min=1e-5, inner_print=false, ϵ_opt=1e-9)
	ibr_opts = IBROptions(Δ_min=1e-9, ibr_iter=100, ordering=[1,2,3,4,5,6,7,8,9,10])

	# Define the game problem
	prob0_alg = GameProblem(N,dt,x0,model,opts_alg,game_obj,game_con)
	prob0_ibr = GameProblem(N,dt,x0,model,opts_ibr,game_obj,game_con)
	prob_alg = deepcopy(prob0_alg)
	prob_ibr = deepcopy(prob0_ibr)

	# Solve the problem
	@elapsed newton_solve!(prob_alg)
	@elapsed ibr_newton_solve!(prob_ibr, ibr_opts=ibr_opts)

	return prob_alg, prob_ibr
end

function display_result(prob::GameProblem; sample::Int=1)
	L = length(prob.stats.res[1:sample:end])
	for t = 1:L
		t_elap = cumsum(prob.stats.t_elap)[1:sample:end][t]
		res = prob.stats.res[1:sample:end][t]
		println("($t_elap, $res)")
	end
	return nothing
end

# Visualize the Results
# Run it twice to avoid including compile times.
plt = plot(layout = (2, 2), xlabel="time (s)", ylabel="residual")
for (i,p) in enumerate([2,4,7,10])
	prob_alg, prob_ibr = ibr_experiment(p)
	plot!(plt[i], cumsum(prob_alg.stats.t_elap), log.(10, prob_alg.stats.res),
		linewidth=3.0, color=:blue, label="Algames")
	plot!(plt[i], cumsum(prob_ibr.stats.t_elap), log.(10, prob_ibr.stats.res),
		linewidth=3.0, color=:orange, label="IBR")

	# plot!(plt, cumsum(prob_alg.stats.t_elap), log.(10, [v.max for v in prob_alg.stats.dyn_vio]))
	# plot!(plt, cumsum(prob_ibr.stats.t_elap), log.(10, [v.max for v in prob_ibr.stats.dyn_vio]))
end
display(plt)


# Script for exporting results to tikz
display_result(prob_alg)
display_result(prob_ibr, sample=10)
