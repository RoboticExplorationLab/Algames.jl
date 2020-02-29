export
    stage_cost,
    cost!,
    cost_expansion,
    cost_gradient,
    cost_hessian,
    cost_gradient!,
    cost_hessian!

"Evaluate the cost at a knot point"
function stage_cost(cost::TO.CostFunction, z::KnotPoint, ind::Vector{Int})
    if TO.is_terminal(z)
        TO.stage_cost(cost, state(z))
    else
        TO.stage_cost(cost, state(z), control(z)[ind])*z.dt
    end
end

"Evaluate the cost for a trajectory (non-allocating)"
@inline function cost!(obj::Objective, Z::Traj, list_ind::Vector{Vector{Int}})
    map!(stage_cost, obj.J, obj.cost, Z, list_ind)
end


function cost_expansion(Q, obj, Z::TO.Traj, pl::Vector{Vector{Int}}, p::Int)
    for i = 1:p
        cost_expansion(Q[i], obj[i], Z, pl[i])
    end
    return nothing
end


"Get Qx, Qu pieces of gradient of cost function, multiplied by dt"
function cost_gradient(cost::TO.CostFunction, z::KnotPoint, ind::Vector{Int})
    Qx, Qu = TO.gradient(cost, state(z), control(z)[ind])
    if TO.is_terminal(z)
        dt_x = 1.0
        dt_u = 0.0
    else
        dt_x = z.dt
        dt_u = z.dt
    end
    return Qx*dt_x, Qu*dt_u
end

"Get Qxx, Quu, Qux pieces of Hessian of cost function, multiplied by dt"
function cost_hessian(cost::TO.CostFunction, z::KnotPoint, ind::Vector{Int})
    Qxx, Quu, Qux = TO.hessian(cost, state(z), control(z)[ind])
    if TO.is_terminal(z)
        dt_x = 1.0
        dt_u = 0.0
    else
        dt_x = z.dt
        dt_u = z.dt
    end
    return Qxx*dt_x, Quu*dt_u, Qux*dt_u
end

function cost_gradient!(E, obj::Objective, Z::Traj, ind::Vector{Int})
    N = length(Z)
    for k in eachindex(Z)
        E.x[k], E.u[k] = cost_gradient(obj[k], Z[k], ind)
    end
end

function cost_hessian!(E, obj::Objective, Z::Traj, ind::Vector{Int})
    N = length(Z)
    for k in eachindex(Z)
        E.xx[k], E.uu[k], E.ux[k] = cost_hessian(obj[k], Z[k], ind)
    end
end

"Expand cost for entire trajectory"
function cost_expansion(E, obj::Objective, Z::Traj, ind::Vector{Int})
    cost_gradient!(E, obj, Z, ind)
    cost_hessian!(E, obj, Z, ind)
end


function TO.LQRObjective(Q::Union{Diagonal{T,S},SMatrix}, R::AbstractArray, Qf::AbstractArray, xf::AbstractVector, uf::AbstractVector, N::Int) where {T,S<:SVector}
    n,m = size(Q,1), size(R,1)
    H = @SMatrix zeros(m,n)
    q = -Q*xf
    r = -R*uf
    c = 0.5*xf'*Q*xf + 0.5*uf'*R*uf
    qf = -Qf*xf
    cf = 0.5*xf'*Qf*xf

    ℓ = QuadraticCost(Q, R, H, q, r, c)
    ℓN = QuadraticCost(Qf, R, H, qf, r, cf)
    Objective(ℓ, ℓN, N)
end
