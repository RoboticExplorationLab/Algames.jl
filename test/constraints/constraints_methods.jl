@testset "Constraint Methods" begin

    # Test Collision Avoidance
    T = Float64
    N = 20
    dt = 0.1
    p = 3
    model = DoubleIntegratorGame(p=p)
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)
    radius = 1.0
    add_collision_avoidance!(game_con, radius)

    pu = model.pu
    @test game_con.state_conlist[1].constraints[1].x1 == pu[1]
    @test game_con.state_conlist[1].constraints[2].x1 == pu[1]
    @test game_con.state_conlist[1].constraints[1].x2 == pu[2]
    @test game_con.state_conlist[1].constraints[2].x2 == pu[3]

    @test game_con.state_conlist[2].constraints[1].x1 == pu[2]
    @test game_con.state_conlist[2].constraints[2].x1 == pu[2]
    @test game_con.state_conlist[2].constraints[1].x2 == pu[1]
    @test game_con.state_conlist[2].constraints[2].x2 == pu[3]

    @test game_con.state_conlist[3].constraints[1].x1 == pu[3]
    @test game_con.state_conlist[3].constraints[2].x1 == pu[3]
    @test game_con.state_conlist[3].constraints[1].x2 == pu[1]
    @test game_con.state_conlist[3].constraints[2].x2 == pu[2]

    # Test Spherical Collision Avoidance
    T = Float64
    N = 20
    dt = 0.1
    p = 3
    d = 3
    model = DoubleIntegratorGame(p=p, d=d)
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)
    radius = 1.0
    add_spherical_collision_avoidance!(game_con, radius)

    pu = model.pu
    @test game_con.state_conlist[1].constraints[1].x1 == pu[1]
    @test game_con.state_conlist[1].constraints[2].x1 == pu[1]
    @test game_con.state_conlist[1].constraints[1].x2 == pu[2]
    @test game_con.state_conlist[1].constraints[2].x2 == pu[3]

    @test game_con.state_conlist[2].constraints[1].x1 == pu[2]
    @test game_con.state_conlist[2].constraints[2].x1 == pu[2]
    @test game_con.state_conlist[2].constraints[1].x2 == pu[1]
    @test game_con.state_conlist[2].constraints[2].x2 == pu[3]

    @test game_con.state_conlist[3].constraints[1].x1 == pu[3]
    @test game_con.state_conlist[3].constraints[2].x1 == pu[3]
    @test game_con.state_conlist[3].constraints[1].x2 == pu[1]
    @test game_con.state_conlist[3].constraints[2].x2 == pu[2]

    # Test Control Bounds
    T = Float64
    N = 20
    dt = 0.1
    p = 3
    model = DoubleIntegratorGame(p=p)
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)
    u_min = -10*ones(model.m)
    u_max =  10*ones(model.m)
    u_min[1] = -Inf
    u_max[1] =  Inf
    add_control_bound!(game_con, u_max, u_min)

    @test game_con.control_conlist.constraints[1].u_min == u_min
    @test game_con.control_conlist.constraints[1].u_max == u_max


    # Test Circle Constraint
    T = Float64
    N = 20
    dt = 0.1
    p = 3
    model = DoubleIntegratorGame(p=p)
    px = model.px
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)
    P = 5
    xc = SVector{P,T}([1.0, 2.0, 3.0, 4.0, 5.0])
    yc = SVector{P,T}([-1.0, -2.0, -3.0, -4.0, -5.0])
    radius = SVector{P,T}([0.1, 0.2, 0.3, 0.4, 0.5])
    add_circle_constraint!(game_con, xc, yc, radius)

    @test game_con.state_conlist[1].constraints[1].xi == px[1][1]
    @test game_con.state_conlist[1].constraints[1].yi == px[1][2]

    @test game_con.state_conlist[2].constraints[1].xi == px[2][1]
    @test game_con.state_conlist[2].constraints[1].yi == px[2][2]

    @test game_con.state_conlist[3].constraints[1].xi == px[3][1]
    @test game_con.state_conlist[3].constraints[1].yi == px[3][2]


    # Test WallConstraint
    T = Float64
    N = 20
    dt = 0.1
    p = 3
    model = DoubleIntegratorGame(p=p)
    px = model.px
    n = model.n
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)
    P = 5
    x = 4
    y = 2
    x1 = SVector{P,T}([ 0.0,  0.0,  1.0,  3.0, -2.0])
    y1 = SVector{P,T}([ 1.0, -1.0,  2.0,  2.0,  0.0])
    x2 = SVector{P,T}([ 1.0,  1.0,  2.0,  2.0,  0.0])
    y2 = SVector{P,T}([ 0.0,  0.0,  1.0,  1.0,  0.0])
    xv = SVector{P,T}([ 1.0,  1.0,  1.0,  1.0,  0.0])./sqrt(2)
    yv = SVector{P,T}([ 1.0, -1.0,  1.0, -1.0,  sqrt(2)])./sqrt(2)

    con = WallConstraint(n,x1,y1,x2,y2,xv,yv,x,y)
    walls = [Wall([x1[j],y1[j]], [x2[j],y2[j]], [xv[j],yv[j]]) for j=1:P]
    add_wall_constraint!(game_con, walls )

    @test game_con.state_conlist[1].constraints[1].x == px[1][1]
    @test game_con.state_conlist[1].constraints[1].y == px[1][2]

    @test game_con.state_conlist[2].constraints[1].x == px[2][1]
    @test game_con.state_conlist[2].constraints[1].y == px[2][2]

    @test game_con.state_conlist[3].constraints[1].x == px[3][1]
    @test game_con.state_conlist[3].constraints[1].y == px[3][2]


    # Test helpers
    T = Float64
    N = 20
    dt = 0.1
    p = 3
    model = DoubleIntegratorGame(p=p)
    px = model.px
    n = model.n
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)

    # Add control constraint
    u_min = -10*ones(model.m)
    u_max =  10*ones(model.m)
    add_control_bound!(game_con, u_max, u_min)

    # Add state constraint
    P = 5
    xc = SVector{P,T}([1.0, 2.0, 3.0, 4.0, 5.0])
    yc = SVector{P,T}([-1.0, -2.0, -3.0, -4.0, -5.0])
    radius = SVector{P,T}([0.1, 0.2, 0.3, 0.4, 0.5])
    add_circle_constraint!(game_con, xc, yc, radius)

    # Set constraint parameters
    opts = Options()
    opts.ρ_0 = 1e-3
    opts.ρ_increase = 1e1
    opts.ρ_max = 1e-1
    opts.λ_max = 1e1
    set_constraint_params!(game_con, opts)

    # Test penaly_update and reset_penalties
    @test game_con.control_conval[1].params.μ0 == 1e-3
    @test game_con.state_conval[1][1].params.μ0 == 1e-3

    @test game_con.control_conval[1].params.ϕ == 1e1
    @test game_con.state_conval[1][1].params.ϕ == 1e1

    @test game_con.control_conval[1].params.μ_max == 1e-1
    @test game_con.state_conval[1][1].params.μ_max == 1e-1

    reset_penalties!(game_con)
    @test game_con.control_conval[1].μ[1] == 1e-3*ones(2model.m)
    @test game_con.state_conval[1][1].μ[1] == 1e-3*ones(P)

    penalty_update!(game_con)
    @test game_con.control_conval[1].μ[1] == 1e-2*ones(2model.m)
    @test game_con.state_conval[1][1].μ[1] == 1e-2*ones(P)

    penalty_update!(game_con)
    @test game_con.control_conval[1].μ[1] == 1e-1*ones(2model.m)
    @test game_con.state_conval[1][1].μ[1] == 1e-1*ones(P)

    penalty_update!(game_con)
    penalty_update!(game_con)
    penalty_update!(game_con)
    penalty_update!(game_con)
    @test game_con.control_conval[1].μ[1] == 1e-1*ones(2model.m)
    @test game_con.state_conval[1][1].μ[1] == 1e-1*ones(P)

    reset_penalties!(game_con)
    @test game_con.control_conval[1].μ[1] == 1e-3*ones(2model.m)
    @test game_con.state_conval[1][1].μ[1] == 1e-3*ones(P)


    # Test dual_update and reset_duals
    @test game_con.control_conval[1].params.λ_max == 1e1
    @test game_con.state_conval[1][1].params.λ_max == 1e1

    reset_penalties!(game_con)
    reset_duals!(game_con)
    @test game_con.control_conval[1].λ[1] == 0.0*ones(2model.m)
    @test game_con.state_conval[1][1].λ[1] == 0.0*ones(P)

    dual_update!(game_con)
    @test game_con.control_conval[1].λ[1] == 0.0*ones(2model.m)
    @test game_con.state_conval[1][1].λ[1] == 0.0*ones(P)

    pdtraj = PrimalDualTraj(probsize, dt)
    init_traj!(pdtraj, f=(rng,args)->ones(args), amplitude=1e2)
    evaluate!(game_con, pdtraj.pr)
    @test game_con.control_conval[1].vals[1] == [90*ones(model.m); -110*ones(model.m)]
    dual_update!(game_con)
    @test game_con.control_conval[1].λ[1] == 1e-3 * [90*ones(model.m); 0*ones(model.m)]
    @test game_con.state_conval[1][1].λ[1] == 1e-3 * max.(0, game_con.state_conval[1][1].vals[1])

    init_traj!(pdtraj, f=(rng,args)->ones(args), amplitude=1e5)
    evaluate!(game_con, pdtraj.pr)
    @test game_con.control_conval[1].vals[1] == [(1e5-10)*ones(model.m); -(1e5+10)*ones(model.m)]
    dual_update!(game_con)
    @test game_con.control_conval[1].λ[1] == [opts.λ_max*ones(model.m); 0*ones(model.m)]

    reset_duals!(game_con)
    @test game_con.control_conval[1].λ[1] == 0.0*ones(2model.m)
    @test game_con.state_conval[1][1].λ[1] == 0.0*ones(P)

end
