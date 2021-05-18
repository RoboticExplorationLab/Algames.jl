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


    # Test Wall3DConstraint
    T = Float64
    n = 6
    P = 5
    x = 4
    y = 2
    z = 1
    X  =1
    X0 = SVector{n,T}([0.0, 0.10, -12.0, 0.10, 12.0, 11.0])
    X1 = SVector{n,T}([1.0, 0.55, -12.0, 0.55, 12.0, 11.0])
    X2 = SVector{n,T}([1.0, 1.25, -12.0, 0.75, 12.0, 11.0])
    x1 = SVector{P,T}([ 0.0,  0.0,  0.0,  0.0,  0.0])
    y1 = SVector{P,T}([ 0.0,  0.0,  0.0,  0.0,  0.0])
    z1 = SVector{P,T}([ 0.0,  0.0,  0.0,  0.0,  0.0])

    x2 = SVector{P,T}([ 1.0,  1.0,  1.0,  1.0,  1.0])
    y2 = SVector{P,T}([ 0.0,  0.0,  0.0,  0.0,  0.0])
    z2 = SVector{P,T}([ 0.0,  0.0,  1.0,  0.0,  0.0])

    x3 = SVector{P,T}([ 1.0,  1.0,  1.0,  0.0,  0.0])
    y3 = SVector{P,T}([ 1.0,  1.0,  1.0,  1.0,  1.0])
    z3 = SVector{P,T}([ 0.0,  1.0,  1.0,  0.0,  1.0])

    xv = SVector{P,T}([ 0.0,  0.0,         -1.0/sqrt(2),  0.0,  0.0])
    yv = SVector{P,T}([ 0.0, -1.0/sqrt(2),  0.0,          0.0, -1.0/sqrt(2)])
    zv = SVector{P,T}([ 1.0,  1.0/sqrt(2),  1.0/sqrt(2),  1.0,  1.0/sqrt(2)])

    con = Wall3DConstraint(n,x1,y1,z1,x2,y2,z2,x3,y3,z3,xv,yv,zv,x,y,z)
    @test norm(TrajectoryOptimization.evaluate(con,X0) - [0.0, -0.1/sqrt(2), -0.1/sqrt(2), 0.0, -0.1/sqrt(2)], 1) < 1e-10
    @test norm(TrajectoryOptimization.evaluate(con,X1) - [1.0,  0.45/sqrt(2),  0.45/sqrt(2), 1.0,  0.45/sqrt(2)], 1) < 1e-10
    @test norm(TrajectoryOptimization.evaluate(con,X2) - [0.0,  0.0,  0.0,  1.0,  -0.25/sqrt(2)], 1) < 1e-10
    @test (@ballocated TrajectoryOptimization.evaluate($con,$X0)) == 0

    ∇c_easy = easy_jacobian(con, X0)
    ∇c = zeros(MMatrix{P,n,T,P*n})
    jacobian!(∇c,con,X0)
    @test norm(∇c - ∇c_easy, 1) < 1e-10
    @test (@ballocated jacobian!($∇c,$con,$X0)) == 0
end
