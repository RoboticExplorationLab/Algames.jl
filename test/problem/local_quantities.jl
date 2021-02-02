@testset "Local quantities" begin

    # Test Dynamics
    T = Float64
    dt = 0.01
    model = DoubleIntegratorGame(p=3, d=2)
    x = rand(SVector{model.n,T})
    u = rand(SVector{model.m,T})
    x_dot = RobotDynamics.dynamics(model, x, u)
    x1 = x + dt*x_dot

    z = KnotPoint(x,u,dt)
    @test typeof(z) <: AbstractKnotPoint
    @test norm(dynamics_residual(model, z, x1), 1) < 1e-3
    @test (@ballocated dynamics_residual($model, $z, $x1)) == 0

    N = 10
    k = 5
    probsize = ProblemSize(N,model)
    pdtraj = PrimalDualTraj(probsize, dt)
    @test (@ballocated dynamics_residual($model, $pdtraj, $k)) == 0

    # Test dynamics Jacobian
    function ∇dynamics_easy(∇f, model, z)
        dt = z.dt
        x = state(z)
        u = control(z)
        function local_dyn_x(x)
            z = KnotPoint(x, u, dt)
            return discrete_dynamics(RK2, model, z)
        end
        function local_dyn_u(u)
            z = KnotPoint(x, u, dt)
            return discrete_dynamics(RK2, model, z)
        end
        ∇f[:,1:n] = ForwardDiff.jacobian(local_dyn_x, x)
        ∇f[:,n .+ (1:m)] = ForwardDiff.jacobian(local_dyn_u, u)
        return ∇f
    end

    T = Float64
    dt = 0.2
    model = UnicycleGame(p=3)
    x = rand(SVector{model.n,T})
    u = rand(SVector{model.m,T})
    z = KnotPoint(x,u,dt)

    n = model.n
    m = model.m
    ∇f1 = zeros(MMatrix{n,n+m,T,n*(n+m)})
    ∇f2 = zeros(MMatrix{n,n+m,T,n*(n+m)})
    ∇f3 = zeros(MMatrix{n,n+m,T,n*(n+m)})
    RobotDynamics.discrete_jacobian!(RK2,∇f1,model,z)
    ∇dynamics_easy(∇f2, model, z)
    ∇dynamics!(∇f3, model, z)
    @test norm(∇f1 - ∇f2, 1) < 1e-10
    @test norm(∇f1 - ∇f3, 1) < 1e-10
    @test (@ballocated ∇dynamics!($∇f3, $model, $z)) == 0
    @test (@ballocated ∇dynamics!($∇f3, $model, $pdtraj, $k)) == 0

end
