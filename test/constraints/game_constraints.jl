@testset "Game Constraints" begin

    # Test GameConstraintValues
    N = 10
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    game_con = GameConstraintValues(probsize)
    @test game_con.probsize.p == model.p
    @test length(game_con.state_conlist) == model.p

    # Test set_constraint_params!
    T = Float64
    u_max =  rand(SVector{model.m,T})
    u_min = -rand(SVector{model.m,T})
    add_control_bound!(game_con, u_max, u_min)
    radius = 1.0
    add_collision_avoidance!(game_con, radius)
    opts = Options{T}()
    opts.ρ_increase = 2.0
    opts.ρ_0 = 3.0
    opts.ρ_max = 4.0
    opts.λ_max = 5.0
    set_constraint_params!(game_con, opts)

    @test game_con.α_dual == opts.α_dual
    @test game_con.αx_dual == opts.αx_dual[1:p]
    @test game_con.active_set_tolerance == opts.active_set_tolerance

    @test game_con.state_conval[model.p][1].params.ϕ == 2.0
    @test game_con.state_conval[model.p][1].params.μ0 == 3.0
    @test game_con.state_conval[model.p][1].params.μ_max == 4.0
    @test game_con.state_conval[model.p][1].params.λ_max == 5.0

    @test game_con.control_conval[1].params.ϕ == 2.0
    @test game_con.control_conval[1].params.μ0 == 3.0
    @test game_con.control_conval[1].params.μ_max == 4.0
    @test game_con.control_conval[1].params.λ_max == 5.0


end
