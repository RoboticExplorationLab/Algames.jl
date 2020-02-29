export
	add_dynamics_constraints!,
	StaticInds

# "```julia
# add_dynamics_constraints!(prob::Problem)
# ```
# Add dynamics constraints to the constraint set"
function add_dynamics_constraints!(dyn_conSet::ConstraintSet{T}, prob::GameProblem{Q}) where {T,Q}
	# n,m = size(prob)
    # Implicit dynamics
    dyn_con = ConstraintVals( DynamicsConstraint{Q}(prob.model, prob.N), 1:prob.N-1 )
    add_constraint!(dyn_conSet, dyn_con, 1)

    # # Initial condition
    # init_con = ConstraintVals( GoalConstraint(prob.x0), 1:1)
    # add_constraint!(conSet, init_con, 1)

    return nothing
end

struct StaticInds{ln,lm,lnm,l2nm,l2n2m}
	s0::SVector{0,Int}
	sn::SVector{ln,Int}
	sm::SVector{lm,Int}
	snm::SVector{lnm,Int}
	s2nm::SVector{l2nm,Int}
	s2n2m::SVector{l2n2m,Int}
end
