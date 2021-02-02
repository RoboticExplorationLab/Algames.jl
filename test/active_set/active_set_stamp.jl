@testset "Constraint Stamp" begin

    # CStamp
    N = 10
    p = 4
    s0 = CStamp(:v, :col, 1, 2, 3)
    @test valid(s0, N, p)
    s0 = CStamp(:v, :col, 1, 1, 3)
    @test !valid(s0, N, p)
    s0 = CStamp(:v, :col, 1, 3, 3)
    @test valid(s0, N, p)
    s0 = CStamp(:v, :col, 1, 5, 3)
    @test !valid(s0, N, p)
    s0 = CStamp(:v, :col, 0, 3, 3)
    @test !valid(s0, N, p)
    s0 = CStamp(:v, :col, 1, 3, 1)
    @test !valid(s0, N, p)
    s0 = CStamp(:v, :col, 3, 1, 3)
    @test !valid(s0, N, p)
    s0 = CStamp(:v, :col, 1, 2, 11)
    @test !valid(s0, N, p)

    s0 = CStamp(:h, :col, 1, 3, 3)
    @test valid(s0, N, p)
    s0 = CStamp(:h, :col, 2, 2, 3)
    @test !valid(s0, N, p)
    s0 = CStamp(:h, :col, 3, 1, 3)
    @test valid(s0, N, p)
    s0 = CStamp(:h, :col, 3, 1, 1)
    @test !valid(s0, N, p)
    s0 = CStamp(:h, :col, 3, 1, 11)
    @test !valid(s0, N, p)
    s0 = CStamp(:h, :col, 5, 1, 10)
    @test !valid(s0, N, p)
    s0 = CStamp(:h, :col, 4, 1, 10)
    @test valid(s0, N, p)

    s1 = CStamp()
    @test s1 != s0
    @test !(s1 == s0)

    s2 = stampify(:v,:col,1,2,10)
    stampify!(s1,:v,:col,1,2,10)
    @test s2 == s1
    @test isequal(s2, s1)
    @test hash(s2) == hash(s1)


end
