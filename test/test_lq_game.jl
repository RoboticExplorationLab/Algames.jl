using Test

atol = 1e-14
# Solving a linear quadratic game with ALGAMES
#  - Should be done in one outer iteration
#  - Should be done in 2 inner iterations (one to take a step the second one to check we are at the optimal point)
#  - Should converge perfectly: optimalty constraint ≈ 0.0
algames_solver = GameProblems.algames_linear_quadratic_solver
reset!(algames_solver, reset_type=:full)
solve!(algames_solver)
iter = algames_solver.stats.iterations
inner_iter = algames_solver.stats.iterations_inner[iter]
@test iter == 1
@test inner_iter == 2
@test isapprox(0., algames_solver.stats.optimality_merit[iter][inner_iter], atol=atol)


# Solving a linear quadratic game with ilQGames
#  - Should be done in 2 inner iterations (one to take a step the second one to check we are at the optimal point)
#  - Should converge perfectly: optimalty constraint ≈ 0.0
ilqgames_solver = GameProblems.ilqgames_linear_quadratic_solver
reset!(ilqgames_solver, reset_type=:full)
ilqgames_solver.opts.line_search_upper_bound=Inf

solve!(ilqgames_solver)
iter = ilqgames_solver.stats.iterations
@test iter == 2
@test isapprox(0., ilqgames_solver.stats.gradient[iter], atol=atol)
