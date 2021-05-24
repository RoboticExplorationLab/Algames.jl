@testset "Statistics" begin

      # Test Statistics
      stats = Statistics()
      T = Float64
      N = 10
      dt = 0.1
      p = 3
      model = UnicycleGame(p=p)
      probsize = ProblemSize(N,model)
      x0 = rand(SVector{probsize.n})

      k = 1
      t_elap = 0.2
      res_norm = 0.3
      Δ_traj = 0.1
      dyn_vio = DynamicsViolation(N)
      con_vio = ControlViolation(N)
      sta_vio = StateViolation(N)
      opt_vio = OptimalityViolation(N)
      record!(stats, t_elap, res_norm, Δ_traj, dyn_vio, con_vio, sta_vio, opt_vio, k)
      @test stats.iter == 1

      game_con = GameConstraintValues(probsize)
      u_max =  0.1*ones(model.m)
      u_min = -0.1*ones(model.m)
      add_control_bound!(game_con, u_max, u_min)
      walls = [Wall([0.,1], [1,0], [1,1]/sqrt(2))]
      add_wall_constraint!(game_con, walls)
      pdtraj = PrimalDualTraj(probsize, dt)
      core = NewtonCore(probsize)

      opts = Options()
      Q = [Diagonal(10*ones(SVector{model.ni[i],T})) for i=1:p] # Quadratic state cost
      R = [Diagonal(0.1*ones(SVector{model.mi[i],T})) for i=1:p] # Quadratic control cost
      # Desrired state
      xf = [SVector{model.ni[1],T}([2,+0.4,0,0]),
            SVector{model.ni[2],T}([2, 0.0,0,0]),
            SVector{model.ni[3],T}([3,-0.4,0,0]),
            ]
      # Desired control
      uf = [zeros(SVector{model.mi[i],T}) for i=1:p]
      # Objectives of the game
      game_obj = GameObjective(Q,R,xf,uf,N,model)
      prob = GameProblem(N,dt,x0,model,opts,game_obj,game_con)

      k = 2
      record!(stats, prob, model, game_con, pdtraj, t_elap, Δ_traj, k)
      @test stats.iter == 2
      @test stats.outer_iter == [1,2]
      @test stats.Δ_traj == [0.1,0.1]
      @test stats.t_elap == [0.2,0.2]
      stats.t_elap

      reset!(stats)
      @test stats.iter == 0
      @test stats.outer_iter == Vector{Int}()

end
