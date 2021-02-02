################################################################################
# Penalty
################################################################################

mutable struct Penalty{T}
	ρ::SVector{1,T}
	ρ_trial::SVector{1,T}
end

function Penalty(ρ::T) where {T}
	return Penalty(SVector{1,T}(ρ), SVector{1,T}(ρ))
end


################################################################################
# GameProblem
################################################################################

mutable struct GameProblem{KN,n,m,T,SVd,SVx}
    probsize::ProblemSize
	model::AbstractGameModel
    core::NewtonCore
	x0::SVx
    game_obj::GameObjective
    game_con::GameConstraintValues
	# colmult::CollisionMultiplier
	pdtraj::PrimalDualTraj{KN,n,m,T,SVd}
	pdtraj_trial::PrimalDualTraj{KN,n,m,T,SVd}
    Δpdtraj::PrimalDualTraj{KN,n,m,T,SVd}
	pen::Penalty{T}
    opts::Options
	stats::Statistics
end

function GameProblem(N::Int, dt::T, x0::SVx, model::AbstractGameModel, opts::Options,
	game_obj::GameObjective, game_con::GameConstraintValues
	) where {T,SVx}

	probsize = ProblemSize(N,model)
	pdtraj = PrimalDualTraj(probsize, dt)
	pdtraj_trial = PrimalDualTraj(probsize, dt)
	Δpdtraj = PrimalDualTraj(probsize, dt)
	core = NewtonCore(probsize)
	# colmult = CollisionMultiplier(probsize)
	stats = Statistics()
	pen = Penalty(opts.ρ_0)
	# Set the constraint parameters according to the options selected for the problem.
	set_constraint_params!(game_con, opts)

	TYPE = (eltype(pdtraj.pr), model.n, model.m, T, eltype(pdtraj.du), typeof(x0))
	return GameProblem{TYPE...}(probsize, model, core, x0, game_obj, game_con, #colmult,
		pdtraj, pdtraj_trial, Δpdtraj, pen, opts, stats)
end
