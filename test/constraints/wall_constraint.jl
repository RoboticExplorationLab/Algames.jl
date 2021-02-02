@testset "Wall Constraint" begin

    # Test WallConstraint
    T = Float64
    n = 4
    P = 5
    x = 4
    y = 2
    X = SVector{n,T}([13.0, 1.0, -12.0, 1.0])
    x1 = SVector{P,T}([ 0.0,  0.0,  1.0,  3.0, -2.0])
    y1 = SVector{P,T}([ 1.0, -1.0,  2.0,  2.0,  0.0])
    x2 = SVector{P,T}([ 1.0,  1.0,  2.0,  2.0,  0.0])
    y2 = SVector{P,T}([ 0.0,  0.0,  1.0,  1.0,  0.0])
    xv = SVector{P,T}([ 1.0,  1.0,  1.0,  1.0,  0.0])./sqrt(2)
    yv = SVector{P,T}([ 1.0, -1.0,  1.0, -1.0,  sqrt(2)])./sqrt(2)

    con = WallConstraint(n,x1,y1,x2,y2,xv,yv,x,y)
    @test norm(TrajectoryOptimization.evaluate(con,X) - [sqrt(2)/2, 0.0, -sqrt(2)/2, 0.0, 0.0], 1) < 1e-10
    @test (@ballocated TrajectoryOptimization.evaluate($con,$X)) == 0

    function easy_jacobian(con,X)
        function local_evaluate(X)
            return TrajectoryOptimization.evaluate(con,X)
        end
        return ForwardDiff.jacobian(local_evaluate, X)
    end
    ∇c_easy = easy_jacobian(con, X)
    ∇c = zeros(MMatrix{P,n,T,P*n})
    jacobian!(∇c,con,X)
    @test norm(∇c - ∇c_easy, 1) < 1e-10
    @test (@ballocated jacobian!($∇c,$con,$X)) == 0

end
