@testset "Constraint Derivatives" begin

    # Test gradient and Hessian of constraints
    T = Float64
    N = 10
    dt = 0.1
    model = UnicycleGame(p=3)
    probsize = ProblemSize(N,model)
    pdtraj = PrimalDualTraj(probsize, dt, f=(rng,args)->ones(args), amplitude=0.1)
    game_con = GameConstraintValues(probsize)

    u_max =  ones(SVector{model.m,T})
    u_min = -ones(SVector{model.m,T})
    add_control_bound!(game_con, u_max, u_min)

    cval = game_con.control_conval[1]
    # Change the duals
    for k = 1:N-2
        cval.λ[k] = k*ones(2model.m)
    end
    # Evaluate the constraint and its Jacobian
    TrajectoryOptimization.evaluate!(cval, pdtraj.pr)
    @test cval.vals[1] == [-0.9*ones(model.m); -1.1*ones(model.m)]
    @test cval.vals[end] == [-0.9*ones(model.m); -1.1*ones(model.m)]
    TrajectoryOptimization.jacobian!(cval, pdtraj.pr)
    @test cval.jac[1] == [I(model.m) ;-I(model.m)]
    @test cval.jac[end] == [I(model.m) ;-I(model.m)]
    TrajectoryOptimization.cost_expansion!(cval)
    Iρ = Diagonal((cval.vals[1] .>= 0) .| (cval.λ[1] .> 0))*cval.μ[1][1]
    @test cval.grad[1] == (cval.jac[1]'*cval.λ[1] + cval.jac[1]'*Iρ*cval.vals[1])
    @test cval.hess[1] == cval.jac[1]'*Iρ*cval.jac[1]
    Iρ = Diagonal((cval.vals[end] .>= 0) .| (cval.λ[end] .> 0))*cval.μ[end][1]
    @test cval.grad[end] == (cval.jac[end]'*cval.λ[end] + cval.jac[end]'*Iρ*cval.vals[end])
    @test cval.hess[end] == cval.jac[end]'*Iρ*cval.jac[end]

end
