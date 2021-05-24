@testset "Violations" begin

    # Test dynamics violation
    T = Float64
    N = 10
    dt = 0.1
    i = 1
    model = UnicycleGame(p=3)
    probsize = ProblemSize(N, model)
    pi = probsize.pz[i]
    pdtraj = PrimalDualTraj(probsize, dt)

    init_traj!(pdtraj, x0=zeros(SVector{model.n,T}), f=zeros)
    dyn_vio = dynamics_violation(model, pdtraj)
    @test dyn_vio.N == pdtraj.probsize.N
    @test dyn_vio.vio == zeros(N-1)
    @test dyn_vio.max == 0.0
    init_traj!(pdtraj, x0=ones(SVector{model.n,T}), f=ones, amplitude=1.0)
    @test dynamics_violation(model, pdtraj).max - maximum(abs.(dynamics_residual(model, pdtraj, 1))) < 1e-10
    @test dynamics_violation(model, pdtraj, i).max - maximum(abs.(dynamics_residual(model, pdtraj, 1)[pi])) < 1e-10

    # Test control violation
    game_con = GameConstraintValues(probsize)
    u_max =  0.1*ones(model.m)
    u_min = -0.1*ones(model.m)
    add_control_bound!(game_con, u_max, u_min)
    init_traj!(pdtraj, x0=zeros(SVector{model.n,T}), f=ones, amplitude=1.0)
    con_vio = control_violation(game_con, pdtraj)
    @test con_vio.N == pdtraj.probsize.N
    @test con_vio.vio == 0.9*ones(N-1)
    @test con_vio.max == 0.9
    con_vio = control_violation(game_con, pdtraj, i)
    @test con_vio.N == pdtraj.probsize.N
    @test con_vio.vio == 0.9*ones(N-1)
    @test con_vio.max == 0.9

    # Test state violation
    game_con = GameConstraintValues(probsize)
    walls = [Wall([0.,1], [1,0], [1,1]/sqrt(2))]
    add_wall_constraint!(game_con, walls)
    init_traj!(pdtraj, x0=zeros(SVector{model.n,T}), f=ones, amplitude=1.0)
    sta_vio = state_violation(game_con, pdtraj)
    @test sta_vio.N == pdtraj.probsize.N
    @test norm(sta_vio.vio - [0; sqrt(2)/2*ones(N-1)...], 1) <= 1e-10
    @test norm(sta_vio.max - sqrt(2)/2,1) < 1e-10
    sta_vio = state_violation(game_con, pdtraj, i)
    @test sta_vio.N == pdtraj.probsize.N
    @test norm(sta_vio.vio - [0; sqrt(2)/2*ones(N-1)...], 1) <= 1e-10
    @test norm(sta_vio.max - sqrt(2)/2,1) < 1e-10

    # Test optimality violation
    core = NewtonCore(probsize)
    opt_vio = optimality_violation(core)
    @test opt_vio.N == core.probsize.N
    @test norm(opt_vio.vio - zeros(N), 1) <= 1e-10
    @test norm(opt_vio.max - 0.0,1) < 1e-10
    opt_vio = optimality_violation(core, i)
    @test opt_vio.N == core.probsize.N
    @test norm(opt_vio.vio - zeros(N), 1) <= 1e-10
    @test norm(opt_vio.max - 0.0,1) < 1e-10

    stamp = stampify(:opt, 2, :u, 2, 5)
    add2sub(core.res_sub[stamp], 1e2*ones(model.mi[2]))
    opt_vio = optimality_violation(core)
    vio = zeros(N)
    vio[5] = 1e2
    @test norm(opt_vio.vio - vio, 1) <= 1e-10
    @test norm(opt_vio.max - 1e2,1) < 1e-10

end
