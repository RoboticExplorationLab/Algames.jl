function constraint_jacobian_residual!(prob::GameProblem, pdtraj::PrimalDualTraj)
    N = prob.probsize.N
    n = prob.probsize.n
    p = prob.probsize.p
    pu = prob.probsize.pu
    game_con = prob.game_con
    jac_sub = prob.core.jac_sub

    stamp = Stamp()
    # State constraints
    for i = 1:p
        for conval in game_con.state_conval[i]
            TrajectoryOptimization.evaluate!(conval, pdtraj.pr)
            TrajectoryOptimization.jacobian!(conval, pdtraj.pr)
            TrajectoryOptimization.cost_expansion!(conval)
            for (j,k) in enumerate(conval.inds)
                stampify!(stamp,  :opt, i, :x, 1, k, :x, 1, k)
                valid(stamp,N,p)  ? add2sub(jac_sub[stamp],  conval.hess[j]) : nothing
            end
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        TrajectoryOptimization.evaluate!(conval, pdtraj.pr)
        TrajectoryOptimization.jacobian!(conval, pdtraj.pr)
        TrajectoryOptimization.cost_expansion!(conval)
        for (j,k) in enumerate(conval.inds)
            for i = 1:p
                stampify!(stamp,  :opt, i, :u, i, k, :u, i, k)
                valid(stamp,N,p)  ? add2sub(jac_sub[stamp],  conval.hess[j][pu[i],pu[i]]) : nothing # Be careful this assumes that the control constraints are decoupled per player
                # To introduce control constraints coupled per player we would need to have a separate set of control constraint for each player.
            end
        end
    end
    return nothing
end


function constraint_residual!(prob::GameProblem, pdtraj::PrimalDualTraj)
    N = prob.probsize.N
    n = prob.probsize.n
    p = prob.probsize.p
    pu = prob.probsize.pu
    game_con = prob.game_con
    res_sub = prob.core.res_sub

    vstamp = VStamp()
    # State constraints
    for i = 1:p
        for conval in game_con.state_conval[i]
            TrajectoryOptimization.evaluate!(conval, pdtraj.pr)
            TrajectoryOptimization.jacobian!(conval, pdtraj.pr)
            TrajectoryOptimization.cost_expansion!(conval)
            for (j,k) in enumerate(conval.inds)
                stampify!(vstamp, :opt, i, :x, 1, k)
                valid(vstamp,N,p) ? add2sub(res_sub[vstamp], conval.grad[j]) : nothing
            end
        end
    end
    # Control constraints
    for conval in game_con.control_conval
        TrajectoryOptimization.evaluate!(conval, pdtraj.pr)
        TrajectoryOptimization.jacobian!(conval, pdtraj.pr)
        TrajectoryOptimization.cost_expansion!(conval)
        for (j,k) in enumerate(conval.inds)
            for i = 1:p
                stampify!(vstamp, :opt, i, :u, i, k)
                valid(vstamp,N,p) ? add2sub(res_sub[vstamp], conval.grad[j][pu[i]]) : nothing # Be careful this assumes that the control constraints are decoupled per player
                # To introduce control constraints coupled per player we would need to have a separate set of control constraint for each player.
            end
        end
    end
    return nothing
end
