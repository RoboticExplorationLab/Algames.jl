@testset "Active Set Methods" begin

    # Test active
    N = 10
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)

    T = Float64
    radius = 1.0
    add_collision_avoidance!(game_con, radius)
    s0 = stampify(:v, :col, 1, 2, 12)
    @test active(game_con, s0)[1] == 0

    game_con.state_conval[1][1].λ[1][1] =  10.0
    game_con.state_conval[1][1].λ[2][1] = -20.0
    game_con.state_conval[1][1].λ[3][1] = -20.0

    game_con.state_conval[1][1].vals[3][1] =  1.0
    game_con.state_conval[1][1].vals[4][1] =  2.0
    game_con.state_conval[1][1].vals[5][1] = -2.0

    Altro.update_active_set!(game_con.state_conval[1][1], Altro.Val(0.0))
    s1 = stampify(:v, :col, 1, 2, 2)
    @test active(game_con, s1)[1] == 1
    s1 = stampify(:v, :col, 1, 2, 3)
    @test active(game_con, s1)[1] == 1
    s1 = stampify(:v, :col, 1, 2, 4)
    @test active(game_con, s1)[1] == 1
    s1 = stampify(:v, :col, 1, 2, 5)
    @test active(game_con, s1)[1] == 1
    s1 = stampify(:v, :col, 1, 2, 6)
    @test active(game_con, s1)[1] == 0


    # Test active_vertical_mask
    N = 10
    dt = 0.1
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    ascore = ActiveSetCore(probsize)
    game_con = GameConstraintValues(probsize)
    radius = 1e-8
    add_collision_avoidance!(game_con, radius)
    n = probsize.n
    m = probsize.m
    S = n*p*(N-1) + m*(N-1) + n*(N-1)
    pdtraj = PrimalDualTraj(probsize, dt)
    # All activated
    active_vertical_mask!(ascore, game_con)
    ascore.vmask == Vector{Int}(1:S+(N-1)*p*(p-1)/2)

    # All deactivated
    Random.seed!(100)
    init_traj!(pdtraj, f=rand, amplitude=1e3)
    update_active_set!(game_con, pdtraj.pr)
    active_vertical_mask!(ascore, game_con)
    ascore.vmask == Vector{Int}(1:S)

    # Test active_horizontal_mask
    N = 10
    dt = 0.1
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    ascore = ActiveSetCore(probsize)
    game_con = GameConstraintValues(probsize)
    radius = 1e-8
    add_collision_avoidance!(game_con, radius)
    n = probsize.n
    m = probsize.m
    S = n*p*(N-1) + m*(N-1) + n*(N-1)
    pdtraj = PrimalDualTraj(probsize, dt)
    # All activated
    active_horizontal_mask!(ascore, game_con)
    ascore.hmask == Vector{Int}(1:S+(N-1)*p*(p-1))

    # All deactivated
    Random.seed!(100)
    init_traj!(pdtraj, f=rand, amplitude=1e3)
    update_active_set!(game_con, pdtraj.pr)
    active_horizontal_mask!(ascore, game_con)
    ascore.hmask == Vector{Int}(1:S)


    # Test Nullspace
    T = Float64
    N = 10
    dt = 0.1
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    x0 = rand(SVector{model.n,T})
    opts = Options()

    Q = [Diagonal(rand(SVector{model.ni[i],T})) for i=1:p]
    R = [Diagonal(rand(SVector{model.mi[i],T})) for i=1:p]
    xf = [i*ones(SVector{model.ni[i],T}) for i=1:p]
    uf = [2i*ones(SVector{model.mi[i],T}) for i=1:p]
    game_obj = GameObjective(Q,R,xf,uf,N,model)

    game_con = GameConstraintValues(probsize)
    radius = 1.0
    add_collision_avoidance!(game_con, radius)
    prob = GameProblem(N, dt, x0, model, opts, game_obj, game_con)

    ascore = ActiveSetCore(probsize)
    residual!(ascore, prob, prob.pdtraj)
    residual_jacobian!(ascore, prob, prob.pdtraj)

    update_nullspace!(ascore, prob, prob.pdtraj)
    @test size(ascore.null.mat) == (probsize.S+(N-1)*p*(p-1), (N-1)*p)
    @test length(ascore.null.vec) == (N-1)*p
    @test length(ascore.null.vec[1]) == probsize.S+(N-1)*p*(p-1)


end
