@testset "Regularizer" begin

    # Test Regularizer
    T = Float64
    reg = Regularizer{T}()
    val = 1e-1
    set!(reg, val)

    @test reg.x == val
    @test reg.u == val
    @test reg.λ == val

    mult!(reg, val)
    @test reg.x == val^2
    @test reg.u == val^2
    @test reg.λ == val^2

end
