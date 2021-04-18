@testset "State Bound Constraint" begin

    # Test ControlBoundConstraint
    T = Float64
    dt = 0.1
    n = 6
    m = 6
    X = SVector{n,T}([13.0, 1.0, -12.0, 1.0, 2.0, 30.0])
    U = SVector{m,T}([13.0, 1.0, -12.0, 1.0, 2.0, 30.0])
    Z = KnotPoint(X,U,dt)
    x_max = [Inf,   Inf,  -11.0, 15.0, 2.0,  30.0]
    x_min = [-Inf, -10.0, -Inf,  1.0, -2.0, -30.0]

    con = StateBoundConstraint(n,x_max=x_max, x_min=x_min)
    @test TrajectoryOptimization.evaluate(con,Z) == [-1.0, -14.0, 0.0, 0.0, -11.0, 0.0, -4.0, -60.0]
    con.inds == [3,4,5,6,8,10,11,12]
    @test (@ballocated TrajectoryOptimization.evaluate($con,$Z)) == 0

    function easy_jacobian(con,Z)
        X = state(Z)
        function local_evaluate(X)
            Z_ = KnotPoint(X, control(Z), Z.dt)
            return TrajectoryOptimization.evaluate(con,Z_)
        end
        return ForwardDiff.jacobian(local_evaluate, state(Z))
    end
    ∇c_easy = easy_jacobian(con, Z)
    ∇c = zeros(MMatrix{8,n,T,8*n})
    jacobian!(∇c,con,Z)
    @test norm(∇c - ∇c_easy, 1) < 1e-10
    @test (@ballocated jacobian!($∇c,$con,$Z)) == 0

end
