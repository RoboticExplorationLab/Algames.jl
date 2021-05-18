@testset "QuadrotorGame" begin

	# Test QuadrotorGame
	p = 3
	model = QuadrotorGame(p=p)
	@test typeof(model) <: QuadrotorGame{model.n, model.m, model.p}
	@test model.n == 12p
	@test model.m == 4p
	@test model.p == p
	@test model.ni == 12*ones(Int,p)
	@test model.mi == 4*ones(Int,p)
	@test model.pu[1] == SVector{4,Int}((1,4,7,10))
	@test model.pu[2] == SVector{4,Int}((2,5,8,11))
	@test model.pu[3] == SVector{4,Int}((3,6,9,12))

	@test model.px[1] == SVector{2,Int}((1,4))
	@test model.px[2] == SVector{2,Int}((2,5))
	@test model.px[3] == SVector{2,Int}((3,6))

	@test model.pz[1] == SVector{12,Int}((1,4,7,10,13,16,19,22,25,28,31,34))
	@test model.pz[2] == SVector{12,Int}((2,5,8,11,14,17,20,23,26,29,32,35))
	@test model.pz[3] == SVector{12,Int}((3,6,9,12,15,18,21,24,27,30,33,36))

	@test size(model) == (12p, 4p, model.pu, p)

	T = Float64
	n = 12
	m = 4

	p = 3
	model = QuadrotorGame(p=p)
	x = rand(SVector{n*p,T})
	u = rand(SVector{m*p,T})
	@test (@ballocated $(Algames.forces)($model, $x, $u, $1)) == 0
	@test (@ballocated $(Algames.moments)($model, $x, $u, $2)) == 0
	@test (@ballocated $(Algames.wrenches)($model, $x, $u, $2)) == 0
	@test (@ballocated RobotDynamics.dynamics($model, $x, $u, $1)) == 0
	for p = 1:4
		model = QuadrotorGame(p=p)
		x = rand(SVector{n*p,T})
		u = rand(SVector{m*p,T})
		@test (@ballocated RobotDynamics.dynamics($model, $x, $u)) == 0
	end

end
