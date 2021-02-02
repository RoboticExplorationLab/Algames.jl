################################################################################
# Dynamics Residual
################################################################################

function dynamics_residual(model::AbstractGameModel, pdtraj::PrimalDualTraj{KN}, k::Int) where {KN}
    zk1 = pdtraj.pr[k+1]::KN
    xk1 = state(zk1)
    zk = pdtraj.pr[k]::KN
    return dynamics_residual(model, zk, xk1)
end

function dynamics_residual(model::AbstractGameModel, zk::AbstractKnotPoint, xk1::SVx) where {SVx}
    return discrete_dynamics(RK2, model, zk) - xk1
end

################################################################################
# Dynamics Jacobian
################################################################################

function ∇dynamics!(∇f::SM, model::AbstractGameModel, pdtraj::PrimalDualTraj{KN}, k::Int) where {KN,SM}
    zk = pdtraj.pr[k]::KN
    return ∇dynamics!(∇f, model, zk)
end

function ∇dynamics!(∇f::SM, model::AbstractGameModel, zk::AbstractKnotPoint) where {SM}
    return RobotDynamics.discrete_jacobian!(RK2, ∇f, model, zk)
end
