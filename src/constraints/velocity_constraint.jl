function add_velocity_bound!(model::AbstractGameModel, game_con::GameConstraintValues, v_max::AbstractVector, v_min::AbstractVector)
	n = model.n
	p = model.p
	@assert length(v_max) == length(v_min) == p
	for i = 1:p
		add_velocity_bound!(model, game_con, i, v_max[i], v_min[i])
	end
	return nothing
end

function add_velocity_bound!(model::AbstractGameModel, game_con::GameConstraintValues, i::Int, v_max::Real, v_min::Real)
	@assert (v_max != Inf) || (v_min != -Inf)
	n = model.n
	p = model.p
	# Compute the bound
	x_max = Inf*ones(n)
	x_min = -Inf*ones(n)
	vel_index_i = velocity_index(model, i)
	x_max[vel_index_i] = v_max
	x_min[vel_index_i] = v_min
	# Add the state bound to all players
	for j = 1:p
		Algames.add_state_bound!(game_con, j, x_max, x_min)
	end
	return nothing
end

function velocity_index(model::UnicycleGame, i::Int)
	@assert i ∈ (1:model.p)
	return model.pz[i][4]
end

function velocity_index(model::BicycleGame, i::Int)
	@assert i ∈ (1:model.p)
	return model.pz[i][3]
end

function velocity_index(model::DoubleIntegratorGame, i::Int)
	error("Velocity Index is not implemented for DoubleIntegratorGame.")
	return nothing
end
