################################################################################
# GameConstraintValues
################################################################################

mutable struct GameConstraintValues{T}
	probsize::ProblemSize
	α_dual::T
	αx_dual::Vector{T}
	active_set_tolerance::T
	state_conlist::Vector{TrajectoryOptimization.AbstractConstraintSet}
	control_conlist::TrajectoryOptimization.AbstractConstraintSet
	state_conval::Vector{Vector{TrajectoryOptimization.AbstractConstraintValues}}
	control_conval::Vector{TrajectoryOptimization.AbstractConstraintValues}
end

function GameConstraintValues(probsize::ProblemSize)
	N = probsize.N
	n = probsize.n
	m = probsize.m
	p = probsize.p
	α_dual = 1.0
	αx_dual = ones(p)
	active_set_tolerance = 0.0
	state_conlist = [ConstraintList(n,m,N) for i=1:p]
	control_conlist = ConstraintList(n,m,N)
	state_conval = [Vector{TrajectoryOptimization.AbstractConstraintValues}() for i=1:p]
	control_conval = Vector{TrajectoryOptimization.AbstractConstraintValues}()
	return GameConstraintValues{typeof(α_dual)}(probsize, α_dual, αx_dual,
		active_set_tolerance, state_conlist,
		control_conlist, state_conval, control_conval)
end

function set_constraint_params!(game_con::GameConstraintValues, opts::Options)
	p = game_con.probsize.p
	game_con.α_dual = opts.α_dual
	game_con.αx_dual = opts.αx_dual[1:p]
	game_con.active_set_tolerance = opts.active_set_tolerance
	for i = 1:game_con.probsize.p
		for conval in game_con.state_conval[i]
			conval.params.ϕ = opts.ρ_increase
			conval.params.μ0 = opts.ρ_0
			conval.params.μ_max = opts.ρ_max
			conval.params.λ_max = opts.λ_max
		end
	end
	for conval in game_con.control_conval
		conval.params.ϕ = opts.ρ_increase
		conval.params.μ0 = opts.ρ_0
		conval.params.μ_max = opts.ρ_max
		conval.params.λ_max = opts.λ_max
	end
	return nothing
end
