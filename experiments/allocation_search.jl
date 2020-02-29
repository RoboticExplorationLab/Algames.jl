include("../src/utils/tests.jl")

k = 1
i = 1
test_allocation(rel_zinds,
    [solver_directgames.constraints.constraints[1],
    k,i,n,m,pl,p,N])

test_allocation(update_g_!, [solver_directgames])
test_allocation(update_g_!, [solver_directgames])

test_allocation(update_H_!, [solver_directgames])
test_allocation(update_H_!, [solver_directgames])

test_allocation(solve!, [solver_directgames])


test_allocation(size, [solver_directgames])
test_allocation(size, [solver_directgames.model])
test_allocation(dual_ascent!, [solver_directgames])
test_allocation(inner_step!, [solver_directgames])
test_allocation(line_search!, [solver_directgames])
test_allocation(primal_dual_update!, [solver_directgames, 1.0])
test_allocation(primal_dual_update!, [solver_directgames, 1.0])
test_allocation(primal_dual_copy_update!, [solver_directgames, 1.0])
test_allocation(primal_dual_copy_update!, [solver_directgames, 1.0])

test_allocation(cost_expansion, [solver_directgames.C,
    solver_directgames.obj, solver_directgames.Z̄,
    solver_directgames.model.pl, solver_directgames.model.p])
test_allocation(cost_expansion, [solver_directgames.C,
    solver_directgames.obj, solver_directgames.Z̄,
    solver_directgames.model.pl, solver_directgames.model.p])

test_allocation(TO.evaluate!, [solver_directgames.constraints,
    solver_directgames.Z])
test_allocation(TO.jacobian!, [solver_directgames.constraints,
    solver_directgames.Z])
test_allocation(TO.evaluate!, [solver_directgames.dyn_constraints,
    solver_directgames.Z])
test_allocation(TO.discrete_jacobian!, [solver_directgames.∇F,
    solver_directgames.model, solver_directgames.Z])
test_allocation(TO.update_active_set!, [solver_directgames.constraints,
    solver_directgames.Z])
# test_allocation(update_g_!, [solver_directgames])
# test_allocation(update_g_!, [solver_directgames])
# test_allocation(update_g_!, [solver_directgames])
# test_allocation(update_g_!, [solver_directgames])
# test_allocation(update_g_!, [solver_directgames])

con_b = BoundaryConstraint(n, B[2],  B[5],  V[1], 1, 2)
con_c = CircleConstraint(n, SVector{1}(1.0), SVector{1}(1.0), SVector{1}(2.0), 1,2)

test_allocation(TO.evaluate!, [solver_directgames.constraints.constraints[1], solver_directgames.Z])
test_allocation(TO.evaluate!, [solver_directgames.constraints.constraints[2], solver_directgames.Z])
test_allocation(TO.evaluate!, [solver_directgames.constraints.constraints[3], solver_directgames.Z])
test_allocation(TO.evaluate!, [solver_directgames.constraints.constraints[4], solver_directgames.Z])
test_allocation(TO.evaluate!, [solver_directgames.constraints.constraints[51], solver_directgames.Z])
test_allocation(TO.evaluate!, [solver_directgames.constraints.constraints[51], solver_directgames.Z])

for i in 1:length(solver_directgames.constraints.constraints)
    @show i
    println(typeof(solver_directgames.constraints.constraints[i].con))
end
