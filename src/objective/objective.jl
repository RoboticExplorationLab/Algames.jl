
################################################################################
# GameObjective
################################################################################

mutable struct GameObjective
	probsize::ProblemSize  # Problem Size
	obj::Vector{Vector{Objective}} # objective of each player [number of players, number of obejctives]
	E::Vector{Vector{Objective}}   # quadratic expansion of the objective of each player [number of players, number of obejctives]
end

function GameObjective(Q::Vector{DQ}, R::Vector{DR}, xf::Vector{SVx}, uf::Vector{SVu},
	N::Int, model::AbstractGameModel) where {DQ,DR,SVx,SVu}
	T = eltype(xf[1])
	n = model.n
	m = model.m
	p = model.p
	probsize = ProblemSize(N,model)
	@assert length(Q) == length(R) == length(xf) == length(uf)

	obj = [Vector{Objective}(undef, 1) for i=1:p]
	E = [Vector{Objective}(undef, 1) for i=1:p]
	for i = 1:p
		Qi = Diagonal(SVector{n,T}(expand_vector(diag(Q[i]), model.pz[i], n)))
		Ri = Diagonal(SVector{m,T}(expand_vector(diag(R[i]), model.pu[i], m)))
		Rf = Diagonal(zeros(SVector{m,T}))
		xfi = SVector{n,T}(expand_vector(xf[i], model.pz[i], n))
		ufi = SVector{m,T}(expand_vector(uf[i], model.pu[i], m))
		cost = LQRCost(Qi,Ri,xfi,ufi,checks=false)
		cost_term = LQRCost(Qi,Rf,xfi,checks=false)
		obj[i][1] = Objective(cost, cost_term, N)
		E[i][1] = Objective([LQRCost(1e-10*ones(n,n)+I, 1e-10*zeros(m,m)+I, zeros(MVector{n,T})) for k=1:N])
	end
	return GameObjective(probsize,obj,E)
end

function expand_vector(v::AbstractVector{T}, inds::AbstractVector{Int}, n::Int) where {T}
	V = zeros(n)
	V[inds] = v
	return V
end

function cost_gradient!(game_obj::GameObjective, pdtraj::PrimalDualTraj)
	p = game_obj.probsize.p
	for i = 1:p
		n_obj = length(game_obj.E[i])
		for j = 1:n_obj
			TrajectoryOptimization.cost_gradient!(game_obj.E[i][j], game_obj.obj[i][j], pdtraj.pr, true)
		end
	end
	return nothing
end

function cost_hessian!(game_obj::GameObjective, pdtraj::PrimalDualTraj)
	p = game_obj.probsize.p
	for i = 1:p
		n_obj = length(game_obj.E[i])
		for j = 1:n_obj
			TrajectoryOptimization.cost_hessian!(game_obj.E[i][j], game_obj.obj[i][j], pdtraj.pr, true, true)
		end
	end
	return nothing
end


function add_collision_cost!(game_obj::GameObjective, radius::AbstractVector{T}, μ::AbstractVector{T}) where {T}
	N = game_obj.probsize.N
	n = game_obj.probsize.n
	m = game_obj.probsize.m
	p = game_obj.probsize.p
	px = game_obj.probsize.px
	@assert p == length(radius) == length(μ)
 	for i = 1:p
		for j ∈ setdiff(1:p,i)
			obj = Objective(CollisionCost{n,m,T,length(px[i])}(μ[i], radius[i], px[i], px[j]), N)
			E = Objective([LQRCost(1e-10*ones(n,n)+I, 1e-10*zeros(m,m)+I, zeros(MVector{n,T})) for k=1:N])
			push!(game_obj.obj[i], obj)
			push!(game_obj.E[i], E)
		end
	end
	return nothing
end




################################################################################
# Collision Cost
################################################################################

mutable struct CollisionCost{n,m,T,ni} <: TrajectoryOptimization.CostFunction
	μ::T
	r::T
	pxi::SVector{ni,Int}
	pxj::SVector{ni,Int}
    terminal::Bool

    function CollisionCost{n,m,T,ni}(μ::T, r::T, pxi::SVector{ni,Int},
		pxj::SVector{ni,Int}; terminal=false) where {n,m,T,ni}
        new{n,m,T,ni}(μ,r,pxi,pxj,terminal)
    end
end

function TrajectoryOptimization.stage_cost(cost::CollisionCost, x::AbstractVector, u::AbstractVector)
    J = TrajectoryOptimization.stage_cost(cost, x)
    return J
end

function TrajectoryOptimization.stage_cost(cost::CollisionCost, x::AbstractVector{T}) where T
	xi = x[cost.pxi]
	xj = x[cost.pxj]
	0.5 * cost.μ * max(0., cost.r - norm(xi-xj))^2
end


function TrajectoryOptimization.gradient!(E::TrajectoryOptimization.QuadraticCostFunction, cost::CollisionCost, x)
	n = length(x)
	ϵ = 1e-10
	ϵ_norm = ϵ*sqrt(n)
	xi = x[cost.pxi]
	xj = x[cost.pxj]
	Δ = (xi-xj)
	Δ_norm = norm(xi-xj)
	E.q .= 0.
	if max(0., cost.r - Δ_norm) > 0.
		g = cost.μ * (cost.r * (ϵ .+ Δ)/(ϵ_norm + Δ_norm) - Δ)
		E.q[cost.pxi] .= - g
	    E.q[cost.pxj] .=   g
	end
    return false
end

function TrajectoryOptimization.gradient!(E::TrajectoryOptimization.QuadraticCostFunction, cost::CollisionCost, x, u)
    TrajectoryOptimization.gradient!(E, cost, x)
    E.r .= 0.
    return false
end

function TrajectoryOptimization.hessian!(E::TrajectoryOptimization.QuadraticCostFunction, cost::CollisionCost, x)
	n = length(x)
	ϵ = 1e-10
	ϵ_norm = ϵ*sqrt(n)
	xi = x[cost.pxi]
	xj = x[cost.pxj]
	Δ = (xi-xj)
	Δ_norm = norm(xi-xj)
	E.Q .= 0.
	if max(0., cost.r - Δ_norm) > 0.
		E.Q[cost.pxi, cost.pxi] .=  cost.μ * (I - cost.r*I/Δ_norm + cost.r*(Δ*Δ')/(Δ_norm^3))
		E.Q[cost.pxi, cost.pxj] .= -cost.μ * (I - cost.r*I/Δ_norm + cost.r*(Δ*Δ')/(Δ_norm^3))
		E.Q[cost.pxj, cost.pxi] .= -cost.μ * (I - cost.r*I/Δ_norm + cost.r*(Δ*Δ')/(Δ_norm^3))
		E.Q[cost.pxj, cost.pxj] .=  cost.μ * (I - cost.r*I/Δ_norm + cost.r*(Δ*Δ')/(Δ_norm^3))
	end
    return true
end

function TrajectoryOptimization.hessian!(E::TrajectoryOptimization.QuadraticCostFunction, cost::CollisionCost, x, u)
    TrajectoryOptimization.hessian!(E, cost, x)
    E.R .= 0
    return true
end

import Base.copy
function Base.copy(c::CollisionCost{n,m,T,ni}) where {n,m,T,ni}
    CollisionCost{n,m,T,ni}(copy(c.μ), copy(c.r), copy(c.pxi), copy(c.pxj), terminal=c.terminal)
end

function TrajectoryOptimization.state_dim(cost::CollisionCost{n,m,T,ni}) where {n,m,T,ni}
	return n
end

function TrajectoryOptimization.control_dim(cost::CollisionCost{n,m,T,ni}) where {n,m,T,ni}
	return m
end
