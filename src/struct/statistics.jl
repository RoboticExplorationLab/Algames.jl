################################################################################
# Statistics
################################################################################

mutable struct Statistics{T}
	iter::Int
	outer_iter::Vector{Int}
	Δ_traj::Vector{T}
	dyn_vio::AbstractVector{DynamicsViolation}
	con_vio::AbstractVector{ControlViolation}
	sta_vio::AbstractVector{StateViolation}
	opt_vio::AbstractVector{OptimalityViolation}
end

function Statistics(;T=Float64)
	iter = 0
	outer_iter = Vector{Int}()
	Δ_traj = zeros(T,0)
	dyn_vio = Vector{DynamicsViolation}()
	con_vio = Vector{ControlViolation}()
	sta_vio = Vector{StateViolation}()
	opt_vio = Vector{OptimalityViolation}()
	return Statistics{T}(iter, outer_iter, Δ_traj, dyn_vio, con_vio, sta_vio, opt_vio)
end

function record!(stats::Statistics, Δ_traj::T, dyn_vio::DynamicsViolation,
	con_vio::ControlViolation, sta_vio::StateViolation, opt_vio::OptimalityViolation, k::Int) where {T}
	stats.iter += 1
	push!(stats.outer_iter, k)
	push!(stats.Δ_traj, Δ_traj)
	push!(stats.dyn_vio, dyn_vio)
	push!(stats.con_vio, con_vio)
	push!(stats.sta_vio, sta_vio)
	push!(stats.opt_vio, opt_vio)
	return nothing
end

function record!(stats::Statistics, core::NewtonCore, model::AbstractGameModel,
	game_con::GameConstraintValues, pdtraj::PrimalDualTraj, Δ_traj::T, k::Int) where {T}

	stats.iter += 1
	push!(stats.outer_iter, k)
	push!(stats.Δ_traj, Δ_traj)
	push!(stats.dyn_vio, dynamics_violation(model, pdtraj))
	push!(stats.con_vio, control_violation(game_con, pdtraj))
	push!(stats.sta_vio, state_violation(game_con, pdtraj))
	push!(stats.opt_vio, optimality_violation(core))
	return nothing
end

function record!(stats::Statistics, core::NewtonCore, model::AbstractGameModel,
	game_con::GameConstraintValues, pdtraj::PrimalDualTraj, Δ_traj::T, k::Int, i::Int) where {T}

	stats.iter += 1
	push!(stats.outer_iter, k)
	push!(stats.Δ_traj, Δ_traj)
	push!(stats.dyn_vio, dynamics_violation(model, pdtraj, i))
	push!(stats.con_vio, control_violation(game_con, pdtraj, i))
	push!(stats.sta_vio, state_violation(game_con, pdtraj, i))
	push!(stats.opt_vio, optimality_violation(core, i))
	return nothing
end

function reset!(stats::Statistics{T}) where {T}
	stats_ = Statistics(T=T)
	for name in fieldnames(Statistics)
		setfield!(stats, name, getfield(stats_, name))
	end
	return nothing
end
