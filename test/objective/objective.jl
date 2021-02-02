@testset "Objective" begin

    # Test GameObjective
    v = [1,2,3]
    inds = [1,3,5]
    n = 5
    V = [1,0,2,0,3]
    @test V == Algames.expand_vector(v,inds,n)

    # Test Game Objective
    T = Float64
    N = 10
    n = 12
    m = 6
    p = 3
    model = UnicycleGame(p=p)
    Q = [Diagonal(rand(SVector{model.ni[i],T})) for i=1:p]
    R = [Diagonal(rand(SVector{model.mi[i],T})) for i=1:p]
    xf = [i*ones(SVector{model.ni[i],T}) for i=1:p]
    uf = [2i*ones(SVector{model.mi[i],T}) for i=1:p]

    game_obj = GameObjective(Q,R,xf,uf,N,model)

    X = [1.0*SVector{model.n,T}([1,2,3,1,2,3,1,2,3,1,2,3]) for k=1:N]
    U = [1.0*SVector{model.m,T}([2,4,6,2,4,6]) for k=1:N-1]
    dt = 0.1
    Dt = dt*ones(N-1)
    traj = Traj(X, U,Dt)

    @test abs(cost(game_obj.obj[1][1], traj)) <= 1e-10
    @test abs(cost(game_obj.obj[2][1], traj)) <= 1e-10
    @test abs(cost(game_obj.obj[3][1], traj)) <= 1e-10

    # Test cost_gradient! and cost_hessian!
    probsize = ProblemSize(N, model)
    pdtraj = PrimalDualTraj(probsize, dt, f=rand, amplitude=1e1)
    for i = 1:p
        for k = 1:N
            game_obj.E[i][1].cost[k].Q += rand(n,n)
            game_obj.E[i][1].cost[k].R += rand(m,m)
            game_obj.E[i][1].cost[k].H += rand(m,n)
            game_obj.E[i][1].cost[k].q += rand(n)
            game_obj.E[i][1].cost[k].r += rand(m)
            game_obj.E[i][1].cost[k].c += rand()
        end
    end
    cost_gradient!(game_obj, pdtraj)
    cost_hessian!(game_obj, pdtraj)

    # Gradient
    # Stage cost
    @test norm(game_obj.E[1][1].cost[1].q - (game_obj.obj[1][1].cost[1].Q*state(pdtraj.pr[1]) + game_obj.obj[1][1].cost[1].q)*dt, 1) < 1e-10
    @test norm(game_obj.E[1][1].cost[1].r - (game_obj.obj[1][1].cost[1].R*control(pdtraj.pr[1]) + game_obj.obj[1][1].cost[1].r)*dt, 1) < 1e-10
    # Terminal cost
    @test norm(game_obj.E[1][1].cost[end].q - (game_obj.obj[1][1].cost[end].Q*state(pdtraj.pr[end]) + game_obj.obj[1][1].cost[end].q), 1) < 1e-10
    @test norm(game_obj.E[1][1].cost[end].r, 1) < 1e-10

    # Hessian
    # Stage cost
    @test norm(game_obj.E[1][1].cost[1].Q - game_obj.obj[1][1].cost[1].Q*dt, 1) < 1e-10
    @test norm(game_obj.E[1][1].cost[1].R - game_obj.obj[1][1].cost[1].R*dt, 1) < 1e-10
    # Terminal cost
    @test norm(game_obj.E[1][1].cost[end].Q - game_obj.obj[1][1].cost[end].Q, 1) < 1e-10
    @test norm(game_obj.E[1][1].cost[end].R, 1) < 1e-10

    # Test add_collision_cost!
    T = Float64
    dt = 0.1
    N = 10
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    n = model.n
    m = model.m
    Q = [Diagonal(rand(SVector{model.ni[i],T})) for i=1:p]
    R = [Diagonal(rand(SVector{model.mi[i],T})) for i=1:p]
    xf = [i*ones(SVector{model.ni[i],T}) for i=1:p]
    uf = [2i*ones(SVector{model.mi[i],T}) for i=1:p]

    game_obj = GameObjective(Q,R,xf,uf,N,model)

    radius = [1.0, 2.0, 3.0]
    μ = [10.0, 20.0, 30.0]
    i = 1
    j = 2
    c = CollisionCost{n,m,T,length(model.px[i])}(μ[i], radius[i], model.px[i], model.px[j])
    add_collision_cost!(game_obj, radius, μ)

    length(game_obj.obj[1]) == 1 + p-1
    length(game_obj.obj[2]) == 1 + p-1
    length(game_obj.obj[3]) == 1 + p-1
    add_collision_cost!(game_obj, radius, μ)
    length(game_obj.obj[1]) == 1 + 2*(p-1)
    length(game_obj.obj[2]) == 1 + 2*(p-1)
    length(game_obj.obj[3]) == 1 + 2*(p-1)

    typeof(game_obj.obj[1][2].cost) <: Vector{<:CollisionCost}
    game_obj.obj[1][2].cost[1].μ == 10.0
    game_obj.obj[2][2].cost[1].μ == 20.0
    game_obj.obj[3][2].cost[1].μ == 30.0

    game_obj.obj[1][2].cost[1].r == 1.0
    game_obj.obj[2][2].cost[1].r == 2.0
    game_obj.obj[3][2].cost[1].r == 3.0



    # Test CollisionCost
    T = Float64
    dt = 0.1
    model = DoubleIntegratorGame(p=2)
    μ = 10.0
    r = 0.2
    pxi = SVector{length(model.px[1]),Int}(model.px[1])
    pxj = SVector{length(model.px[2]),Int}(model.px[2])
    cs = CollisionCost{model.n,model.m,T,length(model.px[1])}(μ,r,pxi,pxj)

    cs0 = copy(cs)
    @test cs0.μ == cs.μ
    @test cs0.r == cs.r
    @test cs0.terminal == cs.terminal
    @test all(cs0.pxi .== cs.pxi)
    @test all(cs0.pxj .== cs.pxj)
    @test state_dim(cs) == model.n
    @test control_dim(cs) == model.m

    @test cs.μ == μ
    @test cs.r == r
    @test cs.pxi == pxi
    @test cs.pxj == pxj
    @test cs.terminal == false
    @test typeof(cs) <: CollisionCost{model.n,model.m,T,length(model.px[1])}
    @test typeof(cs) <: TrajectoryOptimization.CostFunction

    x = @SVector [1.0, 1.1, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    u = rand(SVector{model.m,T})
    z = KnotPoint(x,u,dt)
    @test typeof(z) <: TrajectoryOptimization.AbstractKnotPoint

    # Test Stage Cost
    @test abs(TrajectoryOptimization.stage_cost(cs, x) - 0.05) < 1e-10
    @test abs(TrajectoryOptimization.stage_cost(cs, state(z)) - 0.05) < 1e-10
    @test abs(TrajectoryOptimization.stage_cost(cs, x, u) - 0.05) < 1e-10
    @test abs(TrajectoryOptimization.stage_cost(cs, state(z), control(z)) - 0.05) < 1e-10
    @test abs(TrajectoryOptimization.stage_cost(cs, x) -  0.05) < 1e-10
    @test abs(TrajectoryOptimization.stage_cost(cs, x, u) - 0.05) < 1e-10
    @test abs(TrajectoryOptimization.stage_cost(cs, z) - 0.005) < 1e-10

    @test (@ballocated TrajectoryOptimization.stage_cost($cs, $x)) == 0
    @test (@ballocated TrajectoryOptimization.stage_cost($cs, $x, $u)) == 0
    @test (@ballocated TrajectoryOptimization.stage_cost($cs, $z)) == 0


    # Test Gradient
    function easy_gradient(cs, z)
        function local_eval(x)
            z_ = KnotPoint(x, control(z), z.dt)
            return TrajectoryOptimization.stage_cost(cs,state(z_))
        end
        return ForwardDiff.gradient(local_eval, state(z))
    end

    zr = KnotPoint(rand(SVector{model.n,T}), rand(SVector{model.m,T}), dt)
    cs_active = CollisionCost{model.n,model.m,T,length(model.px[1])}(μ,1e3,pxi,pxj)
    cs_notactive = CollisionCost{model.n,model.m,T,length(model.px[1])}(μ,1e-3,pxi,pxj)
    E = TrajectoryOptimization.QuadraticCost{T}(model.n,model.m)

    TrajectoryOptimization.gradient!(E,cs_active,state(zr))
    grad_x = easy_gradient(cs_active, zr)
    @test norm(E.q - grad_x, 1)/norm(E.q) < 1e-7
    @test E.r == zeros(model.m)

    TrajectoryOptimization.gradient!(E,cs_notactive,state(zr))
    grad_x = easy_gradient(cs_notactive, zr)
    @test norm(E.q - grad_x, 1) < 1e-7
    @test E.r == zeros(model.m)

    @test (@ballocated TrajectoryOptimization.gradient!($E, $cs, $x)) == 0
    @test (@ballocated TrajectoryOptimization.gradient!($E, $cs, $x, $u)) == 0

    # Test Hessian
    function easy_hessian(cs, z)
        function local_eval(x)
            z_ = KnotPoint(x, control(z), z.dt)
            return TrajectoryOptimization.stage_cost(cs,state(z_))
        end
        return ForwardDiff.hessian(local_eval, state(z))
    end

    TrajectoryOptimization.hessian!(E,cs_active,state(zr))
    hess_x = easy_hessian(cs_active, zr)
    E.Q
    hess_x
    @test norm(E.Q - hess_x, 1) < 1e-2
    @test E.r == zeros(model.m)

    TrajectoryOptimization.hessian!(E,cs_notactive,state(zr))
    hess_x = easy_hessian(cs_notactive, zr)
    E.Q
    @test norm(E.Q - hess_x, 1) < 1e-2
    @test E.r == zeros(model.m)

    @test (@ballocated TrajectoryOptimization.hessian!($E, $cs, $x)) == 0
    @test (@ballocated TrajectoryOptimization.hessian!($E, $cs, $x, $u)) == 0

end
