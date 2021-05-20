@testset "Newton Core" begin

    # Test Vertical Indices
    N = 3
    p = 2
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    n = probsize.n
    mi = probsize.mi

    verti_inds = vertical_indices(probsize)
    @test verti_inds[stampify(:opt, 1, :x, 1, 2)] == SVector{n,Int}(1:n)
    @test verti_inds[stampify(:opt, 1, :u, 1, 1)] == SVector{mi[1],Int}(n .+ (1:mi[1]))
    @test verti_inds[stampify(:opt, 1, :x, 1, 3)] == SVector{n,Int}(n + mi[1] .+ (1:n))
    @test verti_inds[stampify(:opt, 1, :u, 1, 2)] == SVector{mi[1],Int}(2n + mi[1] .+ (1:mi[1]))
    @test verti_inds[stampify(:opt, 2, :x, 1, 2)] == SVector{n,Int}(2n+2mi[1] .+ (1:n))

    function test_vertical_indices(probsize::ProblemSize, verti_inds)
        N = probsize.N
        n = probsize.n
        m = probsize.m
        p = probsize.p
        mi = probsize.mi
        ni = probsize.ni
        all_inds = Vector{Int}([])
        for i = 1:p
            for k = 1:N-1
                push!(all_inds, verti_inds[stampify(:opt, i, :x, 1, k+1)]...)
                push!(all_inds, verti_inds[stampify(:opt, i, :u, i, k)]...)
            end
        end
        for k = 1:N-1
            push!(all_inds, verti_inds[stampify(:dyn, 1, :x, 1, k)]...)
        end
        sort!(all_inds)
        valid = true
        valid &= length(all_inds) == p*n*(N-1) + m*(N-1) + n*(N-1)
        valid &= all_inds == [i for i=1:length(all_inds)]
        return valid
    end
    @test test_vertical_indices(probsize, verti_inds)


    # Test Horizontal Indices
    N = 3
    p = 2
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    n = probsize.n
    m = probsize.m
    mi = probsize.mi

    horiz_inds = horizontal_indices(probsize)
    @test horiz_inds[stampify(:x, 1, 2)] == SVector{n,Int}(1:n)
    @test horiz_inds[stampify(:u, 1, 1)] == SVector{mi[1],Int}(n .+ (1:mi[1]))
    @test horiz_inds[stampify(:u, 2, 1)] == SVector{mi[2],Int}(n + mi[1] .+ (1:mi[2]))
    @test horiz_inds[stampify(:λ, 1, 1)] == SVector{n,Int}(n + m .+ (1:n))
    @test horiz_inds[stampify(:λ, 2, 1)] == SVector{n,Int}(2n + m .+ (1:n))
    @test horiz_inds[stampify(:x, 1, 3)] == SVector{n,Int}(3n + m .+ (1:n))

    function test_horizontal_indices(probsize::ProblemSize, horiz_inds)
        N = probsize.N
        n = probsize.n
        m = probsize.m
        p = probsize.p
        mi = probsize.mi
        ni = probsize.ni
        all_inds = Vector{Int}([])
        for k = 1:N-1
            push!(all_inds, horiz_inds[stampify(:x, 1, k+1)]...)
            for i = 1:p
                push!(all_inds, horiz_inds[stampify(:u, i, k)]...)
                push!(all_inds, horiz_inds[stampify(:λ, i, k)]...)
            end
        end
        sort!(all_inds)
        valid = true
        valid &= length(all_inds) == n*(N-1) + m*(N-1) + p*n*(N-1)
        valid &= all_inds == [i for i=1:length(all_inds)]
        return valid
    end
    @test test_horizontal_indices(probsize, horiz_inds)

    # Test dynamics_indices
    N = 10
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N, model)
    dyn = dynamics_indices(probsize)
    n = probsize.n
    mi = probsize.mi
    pu = probsize.pu
    @test dyn[:x][1] == SVector{n,Int}(1:n)
    @test dyn[:u][1] == SVector{mi[1],Int}(n .+ pu[1])
    @test dyn[:u][2] == SVector{mi[2],Int}(n .+ pu[2])
    @test dyn[:u][3] == SVector{mi[3],Int}(n .+ pu[3])

    # Test Newton Core
    N = 3
    p = 2
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    core = NewtonCore(probsize)
    @test typeof(core) <: NewtonCore
    hstamp = stampify(:u, 2, 1)
    vstamp = stampify(:opt, 2, :u, 2, 1)
    stamp = stampify(:opt, 2, :u, 2, 1, :u, 2, 1)
    @test core.horiz_inds[hstamp] == horizontal_idx(core, hstamp)
    @test core.verti_inds[vstamp] == vertical_idx(core, vstamp)
    @test all((core.verti_inds[vstamp], core.horiz_inds[hstamp]) .== idx(core, stamp))

end


@testset "Masks" begin

    # Test Vertical Mask
    N = 6
    p = 5
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    verti_inds = vertical_indices(probsize)
    ni = probsize.ni[1]
    mi = probsize.mi[1]

    msk = []
    for i = 1:p
        m = vertical_mask(probsize, verti_inds, i)
        push!(msk, m)
        @test length(m) == (N - 1) * (2ni + mi)
    end
    @test length(intersect(msk...)) == 0

    # Test horizontal Mask
    N = 6
    p = 5
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N,model)
    horiz_inds = horizontal_indices(probsize)
    ni = probsize.ni[1]
    mi = probsize.mi[1]

    msk = []
    for i = 1:p
        m = horizontal_mask(probsize, horiz_inds, i)
        push!(msk, m)
        @test length(m) == (N - 1) * (2ni + mi)
    end
    @test length(intersect(msk...)) == 0

end
