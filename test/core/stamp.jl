@testset "Stamp" begin
    # Test Stamo and VStamp
    using Test
    N = 10
    p = 3


    vstamp1 = VStamp(:opt, 1, :x, 1, 5)
    @test valid(vstamp1,N,p)

    vstamp1 = VStamp(:opt, 4, :x, 1, 5)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:opt, 2, :u, 2, 10)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:opt, 2, :u, 2, 9)
    @test valid(vstamp1,N,p)

    vstamp1 = VStamp(:opt, 2, :u, 3, 9)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:opt, 2, :z, 1, 1)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:dyn, 1, :u, 2, 1)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:dyn, 1, :x, 1, 1)
    @test valid(vstamp1,N,p)

    vstamp1 = VStamp(:dyn, 2, :x, 1, 1)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:dyn, 2, :x, 2, 1)
    @test !valid(vstamp1,N,p)

    vstamp1 = VStamp(:opt, 2, :u, 2, 4)
    @test valid(vstamp1,N,p)


    hstamp = HStamp()
    @test !valid(hstamp, N, p)

    hstamp1 = HStamp(:u, 2, 3)
    @test valid(hstamp1,N,p)

    hstamp1 = HStamp(:u, 0, 3)
    @test !valid(hstamp1,N,p)

    hstamp1 = HStamp(:u, 1, 100)
    @test !valid(hstamp1,N,p)

    hstamp1 = HStamp(:x, 2, 2)
    @test !valid(hstamp1,N,p)

    hstamp1 = HStamp(:x, 1, 1)
    @test !valid(hstamp1,N,p)

    hstamp1 = HStamp(:x, 1, 2)
    @test valid(hstamp1,N,p)

    hstamp1 = HStamp(:λ, 1, 2)
    @test valid(hstamp1,N,p)

    hstamp1 = HStamp(:λ, 1, 2)
    @test valid(hstamp1,N,p)

    hstamp1 = HStamp(:λ, 0, 2)
    @test !valid(hstamp1,N,p)

    hstamp1 = HStamp(:λ, 2, 4)
    @test valid(hstamp1,N,p)

    hstamp1 = HStamp(:λ, 2, 40)
    @test !valid(hstamp1,N,p)


    stamp1 = Stamp(:opt, 1, :x, 1, 5, :x, 1, 3)
    @test valid(stamp1,N,p)

    stamp1 = Stamp(:opt, 1, :x, 1, 5, :u, 0, 3)
    @test !valid(stamp1,N,p)

    stamp1 = Stamp(:opt, 1, :x, 1, 5, :u, 3, 3)
    @test valid(stamp1,N,p)

    stamp1 = Stamp(:opt, 1, :x, 1, 5, :u, 1, 3)
    @test valid(stamp1,N,p)

    stamp1 = Stamp(:opt, 1, :x, 1, 5, :λ, 3, 3)
    @test !valid(stamp1,N,p)

    stamp1 = Stamp(:opt, 1, :x, 1, 5, :λ, 1, 3)
    @test valid(stamp1,N,p)

    stamp1 = Stamp(:opt, 1, :x, 1, 5, :λ, 1, 3)
    @test valid(stamp1,N,p)

    stamp1 = Stamp(:dyn, 1, :x, 1, 5, :x, 1, 3)
    @test valid(stamp1,N,p)

    stamp1 = Stamp(:dyn, 1, :x, 1, 5, :u, 2, 3)
    @test valid(stamp1,N,p)

    # Test equality
    stamp1 = Stamp(:opt, 1, :x, 1, 5, :λ, 1, 3)
    stamp2 = Stamp(:opt, 1, :x, 1, 5, :λ, 1, 3)
    stamp3 = Stamp(:opt, 2, :x, 1, 5, :λ, 1, 3)
    @test stamp1 == stamp2
    @test !(stamp1 == stamp3)
    @test !(Stamp() == VStamp())
    vstamp1 = VStamp()
    vstamp2 = VStamp()
    vstamp3 = VStamp(:opt, 2, :x, 1, 3)
    @test vstamp1 == vstamp2
    @test !(vstamp1 == vstamp3)
    hstamp1 = HStamp()
    hstamp2 = HStamp()
    hstamp3 = HStamp(:x, 1, 3)
    @test hstamp1 == hstamp2
    @test !(hstamp1 == hstamp3)

    # Test Stampify
    stamp0 = stampify(:opt, 1, :x_1, 1, 6, :λ, 1, 3)
    stamp1 = stampify(:opt, 1, :x, 1, 5, :λ, 1, 3)
    stamp2 = stampify(:opt, 1, :x1, 1, 4, :λ, 1, 3)
    @test stamp0 == stamp1
    @test stamp1 == stamp2

    stampify!(stamp2, :opt, 2, :x1, 1, 4, :λ, 1, 3)
    @test !(stamp1 == stamp2)

    vstamp0 = VStamp()
    vstamp1 = stampify(:opt, 1, :x, 1, 5)
    stampify!(vstamp0, :opt, 1, :x, 1, 5)
    @test vstamp0 == vstamp1
    @test isequal(vstamp0, vstamp1)

    hstamp0 = HStamp()
    hstamp1 = stampify(:x, 1, 5)
    stampify!(hstamp0, :x, 1, 5)
    @test hstamp0 == hstamp1
    @test isequal(hstamp0, hstamp1)

    @test !isequal(vstamp0, hstamp1)

end
