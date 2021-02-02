@testset "Control Bound Constraint" begin

    # Test ControlBoundConstraint
    T = Float64
    dt = 0.1
    n = 4
    m = 6
    X = SVector{n,T}([13.0, 1.0, -12.0, 1.0])
    U = SVector{m,T}([13.0, 1.0, -12.0, 1.0, 2.0, 30.0])
    Z = KnotPoint(X,U,dt)
    u_max = [Inf,   Inf,  -11.0, 15.0, 2.0,  30.0]
    u_min = [-Inf, -10.0, -Inf,  1.0, -2.0, -30.0]

    con = ControlBoundConstraint(m,u_max=u_max, u_min=u_min)
    @test TrajectoryOptimization.evaluate(con,Z) == [-1.0, -14.0, 0.0, 0.0, -11.0, 0.0, -4.0, -60.0]
    con.inds == [3,4,5,6,8,10,11,12]
    @test (@ballocated TrajectoryOptimization.evaluate($con,$Z)) == 0

    function easy_jacobian(con,Z)
        U = control(Z)
        function local_evaluate(U)
            Z_ = KnotPoint(state(Z), U, Z.dt)
            return TrajectoryOptimization.evaluate(con,Z_)
        end
        return ForwardDiff.jacobian(local_evaluate, control(Z))
    end
    ∇c_easy = easy_jacobian(con, Z)
    ∇c = zeros(MMatrix{8,m,T,8*m})
    jacobian!(∇c,con,Z)
    @test norm(∇c - ∇c_easy, 1) < 1e-10
    @test (@ballocated jacobian!($∇c,$con,$Z)) == 0

end
