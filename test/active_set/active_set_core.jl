@testset "Active Set Core" begin

    # ActiveSetCore
    N = 10
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    ascore = ActiveSetCore(probsize)
    @test typeof(ascore) <: ActiveSetCore

end
