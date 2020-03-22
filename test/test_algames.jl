using Test

prob = copy(GameProblems.algames_ramp_merging)
Z0 = copy(prob.Z)
X0 = deepcopy(TO.states(Z0))
U0 = deepcopy(TO.controls(Z0))

# Solve with iLQR
algames = DirectGamesSolver(prob)
solve!(algames)
@test converged(algames)
