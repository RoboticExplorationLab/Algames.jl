@testset "Cylinder Constraint" begin

    # Test CylinderConstraint
    T = Float64
    n = 4
    P = 5
    x = 2
    y = 3
    z = 4
    X0 = SVector{n,T}([13.0, 1.0, 1.0, 2.0])
    p1 = SVector{P,T}([ 1.0,  1.0,  1.0,  0.0,  1.0])
    p2 = SVector{P,T}([ 0.0,  0.0,  0.0,  1.0,  0.0])
    p3 = SVector{P,T}([ 1.0,  1.0,  3.0,  1.0,  2.0])
    v = SVector{P,Symbol}([:z, :z, :z, :x, :y])
    l = SVector{P,T}([ 5.0,  2.0,  0.5,  2.0,  10.0])
    r = SVector{P,T}([ 3.0,  2.0,  7.0,  3.0,  1.0])

    con = CylinderConstraint(n,p1,p2,p3,v,l,r,x,y,z)
    TrajectoryOptimization.evaluate(con,X0)


    @test norm(TrajectoryOptimization.evaluate(con,X0) - [8.0, 3.0, 0.0, 8.0, 1.0], 1) < 1e-10
    @test (@ballocated TrajectoryOptimization.evaluate($con,$X0)) == 0

    function easy_jacobian(con,X)
        function local_evaluate(X)
            return TrajectoryOptimization.evaluate(con,X)
        end
        return ForwardDiff.jacobian(local_evaluate, X)
    end
    ∇c_easy = easy_jacobian(con, X0)
    ∇c = zeros(MMatrix{P,n,T,P*n})
    jacobian!(∇c,con,X0)
    @test norm(∇c - ∇c_easy, 1) < 1e-10
    @test (@ballocated jacobian!($∇c,$con,$X0)) == 0

end
