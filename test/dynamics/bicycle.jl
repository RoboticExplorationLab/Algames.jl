@testset "BicycleGame" begin
    # Test BicycleGame
    p = 3
    model = BicycleGame(p=p,lr=1.0,lf=2.0)
    @test model.lr == 1.0
    @test model.lf == 2.0
    @test typeof(model) <: BicycleGame{model.n, model.m, model.p, model.lf, model.lr}
    @test model.n == 4*p
    @test model.m == 2p
    @test model.p == p
    @test model.ni == 4*ones(Int,p)
    @test model.mi == 2*ones(Int,p)
    @test model.pu[1] == SVector{2,Int}((1,4))
    @test model.pu[2] == SVector{2,Int}((2,5))
    @test model.pu[3] == SVector{2,Int}((3,6))

    @test model.px[1] == SVector{2,Int}((1,4))
    @test model.px[2] == SVector{2,Int}((2,5))
    @test model.px[3] == SVector{2,Int}((3,6))

    @test model.pz[1] == SVector{4,Int}((1,4,7,10))
    @test model.pz[2] == SVector{4,Int}((2,5,8,11))
    @test model.pz[3] == SVector{4,Int}((3,6,9,12))

    @test size(model) == (4*p, 2p, model.pu, p)

    T = Float64
    x = ones(SVector{model.n,T})
    u = 100*ones(SVector{model.m,T})
    @test (@ballocated dynamics($model, $x, $u)) == 0
    # @btime dynamics($model, $x, $u)


end
