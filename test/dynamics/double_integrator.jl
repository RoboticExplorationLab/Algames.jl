@testset "Double Integrator" begin
    # Test DoubleIntegratorGame
    p = 2
    d = 3
    model = DoubleIntegratorGame(p=p, d=d)
    @test model.n == 2d*p
    @test model.m == d*p
    @test model.p == p
    @test model.ni == 2d*ones(Int,p)
    @test model.mi == d*ones(Int,p)
    @test model.pu[1] == SVector{d,Int}((1,3,5))
    @test model.pu[2] == SVector{d,Int}((2,4,6))

    @test model.px[1] == SVector{2,Int}((1,3))
    @test model.px[2] == SVector{2,Int}((2,4))

    @test model.pz[1] == SVector{2d,Int}((1,3,5,7,9,11))
    @test model.pz[2] == SVector{2d,Int}((2,4,6,8,10,12))

    @test size(model) == (2d*p, d*p, model.pu, p)

    T = Float64
    x = ones(SVector{model.n,T})
    u = 100*ones(SVector{model.m,T})
    @test (@ballocated dynamics($model, $x, $u)) == 0
    # @btime dynamics($model, $x, $u)
end
