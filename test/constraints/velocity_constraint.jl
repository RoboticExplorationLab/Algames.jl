@testset "Velocity Constraint" begin

    # UnicycleGame
    T = Float64
    N = 10
    dt = 0.1
    model = UnicycleGame(p=3)
    probsize = ProblemSize(N,model)
    pdtraj = PrimalDualTraj(probsize, dt, f=ones, amplitude=0.1)
    game_con = GameConstraintValues(probsize)

    v_max =  ones(model.p)
    v_min = -ones(model.p)
    add_velocity_bound!(model, game_con, v_max, v_min)

    cval = game_con.state_conval[1]
    @test length(game_con.state_conval) == model.p
    @test length(cval) == model.p

    # BicycleGame
    model = BicycleGame(p=3)
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)

    v_max =  [1,+Inf,+Inf]
    v_min =  [1,-1,-Inf]
    add_velocity_bound!(model, game_con, v_max, v_min)

    cval = game_con.state_conval[1]
    @test length(game_con.state_conval) == model.p
    @test length(cval) == model.p-1

    # Empty constraints
    game_con = GameConstraintValues(probsize)

    v_max =  Inf*ones(model.p)
    v_min = -Inf*ones(model.p)
    add_velocity_bound!(model, game_con, v_max, v_min)

    cval = game_con.state_conval[1]
    @test length(game_con.state_conval) == model.p
    @test length(cval) == 0

end
