################################################################################
# Statistics
################################################################################

mutable struct Statistics
	iter::Int
	outer_iter::Vector{Int}
	dyn_vio::AbstractVector{DynamicsViolation}
	con_vio::AbstractVector{ControlViolation}
	sta_vio::AbstractVector{StateViolation}
	opt_vio::AbstractVector{OptimalityViolation}
end

function Statistics()
	iter = 0
	outer_iter = Vector{Int}()
	dyn_vio = Vector{DynamicsViolation}()
	con_vio = Vector{ControlViolation}()
	sta_vio = Vector{StateViolation}()
	opt_vio = Vector{OptimalityViolation}()
	return Statistics(iter, outer_iter, dyn_vio, con_vio, sta_vio, opt_vio)
end

function record!(stats::Statistics, dyn_vio::DynamicsViolation,
	con_vio::ControlViolation, sta_vio::StateViolation, opt_vio::OptimalityViolation, k::Int)
	stats.iter += 1
	push!(stats.outer_iter, k)
	push!(stats.dyn_vio, dyn_vio)
	push!(stats.con_vio, con_vio)
	push!(stats.sta_vio, sta_vio)
	push!(stats.opt_vio, opt_vio)
	return nothing
end

function record!(stats::Statistics, core::NewtonCore, model::AbstractGameModel,
	game_con::GameConstraintValues, pdtraj::PrimalDualTraj, k::Int)

	stats.iter += 1
	push!(stats.outer_iter, k)
	push!(stats.dyn_vio, dynamics_violation(model, pdtraj))
	push!(stats.con_vio, control_violation(game_con, pdtraj))
	push!(stats.sta_vio, state_violation(game_con, pdtraj))
	push!(stats.opt_vio, optimality_violation(core))
	return nothing
end

function reset!(stats::Statistics)
	stats_ = Statistics()
	for name in fieldnames(Statistics)
		setfield!(stats, name, getfield(stats_, name))
	end
	return nothing
end
