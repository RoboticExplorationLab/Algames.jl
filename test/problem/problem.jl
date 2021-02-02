@testset "Problem" begin

    # Test Penalty
    T = Float64
    pen0 = Penalty{T}(SVector{1,T}(1.0), SVector{1,T}(1.0))
    pen1 = Penalty(1.0)
    @test pen0.ρ == pen1.ρ
    @test pen0.ρ_trial == pen1.ρ_trial


    # Test GameProblem
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
    prob = GameProblem(N, dt, x0, model, opts, game_obj, game_con)
    @test typeof(prob) <: GameProblem

end
