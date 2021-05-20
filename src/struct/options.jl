################################################################################
# Options
################################################################################

@with_kw mutable struct Options{T}
    # Options
	"Gauss-newton convergence criterion: tolerance."
	θ::T=1e-2

	"Initialization function of the primal dual vector."
	f_init::Function=rand

	"Initialization amplitude of the primal dual vector."
	amplitude_init::T=1e-8

	"Shift of the trajectory for the initial guess (useful for MPC, ususally, shift=1 in this case)."
	shift::Int=2^10

	# Regularization
	"Regularization of the residual and residual Jacobian."
	regularize::Bool=true

	"Current Jacobian regularization for each primal-dual variables."
	reg::Regularizer{T}=Regularizer()

	"Initial Jacobian regularization."
	reg_0::T=1e-3

	# Line search
	"Initial step size."
	α_0::T=1.0

	"Line search increase step."
	α_increase::T=1.2

	"Line search decrease step."
	α_decrease::T=0.5

	"Expected residual improvement."
	β::T=0.01

	"Number of line search iterations."
	ls_iter::Int=25

	"Minimum progress allowed before breaking out of inner loop."
	Δ_min::T=1e-9

	# Augmented Lagrangian Penalty
	"Initial augmented Lagrangian penalty."
	ρ_0::T=1.0

	"Fixed augmented Lagrangian penalty on the residual used for the line search trials."
	ρ_trial::T=1.0

	"Penalty increase step."
	ρ_increase::T=10.0

	"Maximum augmented Lagrangian penalty."
	ρ_max::T=1e7

	"Maximum Lagrange multiplier."
	λ_max::T=1e7

	"Control Dual ascent step size for all players."
	α_dual::T=1e0

	"State Dual ascent step size for each player."
	αx_dual::Vector{T}=ones(10)

	"Active set tolerance."
	active_set_tolerance::T=1e-4

	# Constraint satisfaction criterion
	"Dynamics constraint satisfaction criterion."
	ϵ_dyn::T=1e-3

	"State constraint satisfaction criterion."
	ϵ_sta::T=1e-3

	"Control constraint satisfaction criterion."
	ϵ_con::T=1e-3

	"Optimality constraint satisfaction criterion."
	ϵ_opt::T=1e-3

	# Augmented Lagrangian iterations.
	"Outer loop iterations."
	outer_iter::Int=7

	"Inner loop iterations."
	inner_iter::Int=20

	# Problem Scaling
	"Objective function scaling."
	γ::T=1e0

	# MPC
	"Number of time steps simulated with MPC"
	mpc_horizon::Int=20

	"Rate of the upsampling used for simulaion and feedback control in the MPC"
	upsampling::Int=2

	# Printing
	"Displaying the inner step results."
	inner_print::Bool=true

	"Displaying the outer step results."
	outer_print::Bool=true

	"Solver random seed."
	seed::Int=100

	"Reseting the duals."
	dual_reset::Bool=true
end


################################################################################
# Iterative Best Response Solver Options
################################################################################

@with_kw mutable struct IBROptions{T}
    # Options
	"Number of iterative best response iterations."
	ibr_iter::Int=100

	"Ordering of players, in the iterative best response scheme."
	ordering::Vector{Int}=Vector(1:100)

	"Minimum progress allowed before breaking out of the iterative best response loop."
	Δ_min::T=1e-9

	"Live plotting of the trajectories obtained through the iterative best response scheme."
	live_plotting::Bool=false
end
