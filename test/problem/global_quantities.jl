@testset "Global Quantities" begin

    # Test Residual
    T = Float64
    N = 10
    dt = 0.1
    p = 3
    i = 2
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
    prob = GameProblem(N, dt, x0, model, opts, game_obj, game_con)
    @test typeof(prob) <: GameProblem
    residual!(prob)
    ibr_residual!(prob, i)

    # Test Residual Jacobian
    residual_jacobian!(prob)
    ibr_residual_jacobian!(prob, i)

    # Test scn
    @test scn(1234.0) == " 1.2e+3"
    @test scn(-1234.0) == "-1.2e+3"
    @test scn(-0.1234) == "-1.2e-1"
    @test scn( 0.1234) == " 1.2e-1"
    @test scn(0) == " 0.0e+0"
    @test scn(-0) == " 0.0e+0"
    @test scn(0, digits=3) == " 0.000e+0"
    @test scn(1234, digits=3) == " 1.234e+3"
    @test scn(1234, digits=0) == " 1e+3"
    @test_throws AssertionError scn(1234, digits=-1) == " 1e+3"

end
