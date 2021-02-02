@testset "Problem Size" begin

    # Test ProblemSize
    model_dbl = DoubleIntegratorGame(p=3,d=2)
    model_uni = UnicycleGame(p=3)
    N = 10
    probsize_uni = ProblemSize(N,model_uni)
    probsize_dbl = ProblemSize(N,model_dbl)
    @test probsize_uni == probsize_dbl

end
