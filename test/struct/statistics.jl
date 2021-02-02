@testset "Statistics" begin

    # Test Statistics
    stats = Statistics()
    N = 10
    dt = 0.1
    model = UnicycleGame(p=3)
    probsize = ProblemSize(N,model)

    k = 1
    dyn_vio = DynamicsViolation(N)
    con_vio = ControlViolation(N)
    sta_vio = StateViolation(N)
    opt_vio = OptimalityViolation(N)
    record!(stats, dyn_vio, con_vio, sta_vio, opt_vio, k)
    @test stats.iter == 1

    game_con = GameConstraintValues(probsize)
    u_max =  0.1*ones(model.m)
    u_min = -0.1*ones(model.m)
    add_control_bound!(game_con, u_max, u_min)
    walls = [Wall([0.,1], [1,0], [1,1]/sqrt(2))]
    add_wall_constraint!(game_con, walls)
    pdtraj = PrimalDualTraj(probsize, dt)
    core = NewtonCore(probsize)
    k = 2
    record!(stats, core, model, game_con, pdtraj, k)
    @test stats.iter == 2
    @test stats.outer_iter == [1,2]

    reset!(stats)
    @test stats.iter == 0
    @test stats.outer_iter == Vector{Int}()

end
