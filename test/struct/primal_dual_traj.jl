@testset "Primal Dual Traj" begin

    # Test PrimalDualTraj
    N = 10
    dt = 0.2
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N, model)
    pdtraj = PrimalDualTraj(probsize, dt)
    @test size(pdtraj.pr)[1] == N
    @test size(state(pdtraj.pr.data[1]))[1] == model.n
    @test size(control(pdtraj.pr.data[1]))[1] == model.m
    @test pdtraj.pr.data[1].dt == dt
    @test size(pdtraj.du)[1] == p
    @test size(pdtraj.du[1])[1] == N-1
    @test size(pdtraj.du[1][1])[1] == model.n

    # Test init_traj!
    T = Float64
    x0 = rand(SVector{model.n,T})
    init_traj!(pdtraj, x0=x0, f=ones, amplitude=10.0)
    @test state(pdtraj.pr[1]) == x0
    @test state(pdtraj.pr[2]) == 10*ones(model.n)
    @test control(pdtraj.pr[1]) == 10*ones(model.m)
    @test pdtraj.du[1][1] == 10*ones(model.n)
    @test pdtraj.du[1][1] == 10*ones(model.n)

    # Test set_traj!
    n = model.n
    m = model.m
    p = model.p
    core = NewtonCore(probsize)
    x0 = rand(SVector{model.n,T})
    Δpdtraj = PrimalDualTraj(probsize, dt)
    Δtraj = ones(n*(N-1)+m*(N-1)+n*p*(N-1))
    init_traj!(Δpdtraj, x0=x0, f=ones, amplitude=10.0)
    set_traj!(core, Δpdtraj, Δtraj)

    @test state(Δpdtraj.pr[1]) == x0
    @test state(Δpdtraj.pr[2]) == ones(n)
    @test state(Δpdtraj.pr[end]) == ones(n)
    @test control(Δpdtraj.pr[1]) == ones(m)
    @test control(Δpdtraj.pr[end-1]) == ones(m)
    @test Δpdtraj.du[1][1] == ones(n)
    @test Δpdtraj.du[1][end] == ones(n)
    @test Δpdtraj.du[end][end] == ones(n)

    Δtraj = rand(n*(N-1)+m*(N-1)+n*p*(N-1))
    init_traj!(Δpdtraj, x0=x0, f=ones, amplitude=10.0)
    set_traj!(core, Δpdtraj, Δtraj)

    @test state(Δpdtraj.pr[1]) == x0
    @test state(Δpdtraj.pr[2]) == Δtraj[horizontal_idx(core, stampify(:x, 1, 2))]
    @test state(Δpdtraj.pr[end]) == Δtraj[horizontal_idx(core, stampify(:x, 1, N))]
    @test control(Δpdtraj.pr[1])[model.pu[1]] == Δtraj[horizontal_idx(core, stampify(:u, 1, 1))]
    @test control(Δpdtraj.pr[1])[model.pu[2]] == Δtraj[horizontal_idx(core, stampify(:u, 2, 1))]
    @test control(Δpdtraj.pr[2])[model.pu[1]] == Δtraj[horizontal_idx(core, stampify(:u, 1, 2))]
    @test control(Δpdtraj.pr[2])[model.pu[2]] == Δtraj[horizontal_idx(core, stampify(:u, 2, 2))]
    @test control(Δpdtraj.pr[2])[model.pu[3]] == Δtraj[horizontal_idx(core, stampify(:u, 3, 2))]
    @test control(Δpdtraj.pr[end-1])[model.pu[1]] == Δtraj[horizontal_idx(core, stampify(:u, 1, N-1))]
    @test Δpdtraj.du[1][1] == Δtraj[horizontal_idx(core, stampify(:λ, 1, 1))]
    @test Δpdtraj.du[1][end] == Δtraj[horizontal_idx(core, stampify(:λ, 1, N-1))]
    @test Δpdtraj.du[end][end] == Δtraj[horizontal_idx(core, stampify(:λ, 3, N-1))]

    # Test get_traj!
    N = 10
    dt = 0.1
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(n,model)
    core = NewtonCore(probsize)
    Δpdtraj = PrimalDualTraj(probsize, dt)
    Δtraj0 = rand(probsize.S)
    Δtraj1 = deepcopy(Δtraj0)
    Δtraj2 = deepcopy(Δtraj0)

    get_traj!(core, Δpdtraj, Δtraj0)
    Δtraj0
    Δtraj1

    @test !(Δtraj0 == Δtraj1)
    set_traj!(core, Δpdtraj, Δtraj1)
    get_traj!(core, Δpdtraj, Δtraj0)
    @test Δtraj0 == Δtraj1

    # Test update_traj!
    n = model.n
    m = model.m
    p = model.p
    x0 = rand(SVector{model.n,T})
    target = PrimalDualTraj(probsize, dt)
    source = PrimalDualTraj(probsize, dt)
    Δ = PrimalDualTraj(probsize, dt)
    init_traj!(target, x0=x0, f=ones, amplitude=0.0)
    init_traj!(source, x0=x0, f=ones, amplitude=10.0)
    init_traj!(Δ, x0=x0, f=ones, amplitude=100.0)
    α = 0.5
    update_traj!(target, source, α, Δ)

    @test state(target.pr[1]) == x0
    @test state(target.pr[2]) == 60*ones(n)
    @test state(target.pr[end]) == 60*ones(n)
    @test control(target.pr[1]) == 60*ones(m)
    @test control(target.pr[end-1]) == 60*ones(m)
    @test target.du[1][1] == 60*ones(n)
    @test target.du[1][end] == 60*ones(n)
    @test target.du[end][end] == 60*ones(n)

    # Test Δ_step
    T = Float64
    N = 10
    dt = 0.2
    p = 3
    model = UnicycleGame(p=p)
    probsize = ProblemSize(N, model)
    pdtraj = PrimalDualTraj(probsize, dt)
    n = model.n
    m = model.m
    x0 = 1e3*ones(SVector{model.n,T})
    Δ = PrimalDualTraj(probsize, dt)
    init_traj!(pdtraj, x0=x0, f=ones, amplitude=10.0)
    α = 0.5
    @test Δ_step(pdtraj, α) == 10.0*α

end
