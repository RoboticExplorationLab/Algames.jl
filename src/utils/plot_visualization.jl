export
    visualize_trajectory_car,
    visualize_control,
    visualize_state,
    visualize_collision_avoidance,
    visualize_circle_collision,
    visualize_boundary_collision,
    visualize_dynamics,
    visualize_optimality_merit,
    visualize_H_cond,
    visualize_α,
    visualize_cmax,
    circle_shape,
    rectangle_shape


function visualize_trajectory_car(solver::TO.ALTROSolver{T}; title::String="",
    plt::Plots.Plot=plot(), save_figure::Bool=false) where T

    conSet = solver.solver_al.solver_uncon.obj.constraints
    x0 = solver.solver_al.solver_uncon.x0
    xf = solver.solver_al.solver_uncon.xf
    visualize_trajectory_car(solver, x0, xf, conSet; title=title, plt=plt, save_figure=save_figure)
    return nothing
end

function visualize_trajectory_car(solver::TO.AbstractSolver{T}; title::String="",
    plt::Plots.Plot=plot(), save_figure::Bool=false) where T

    conSet = solver.constraints
    x0 = solver.x0
    xf = solver.xf
    visualize_trajectory_car(solver, x0, xf, conSet; title=title, plt=plt, save_figure=save_figure)
    return nothing
end

function visualize_trajectory_car(solver::TO.AbstractSolver{T}, x0::AbstractArray,
    xf::SVector{n1,T}, conSet::TO.ConstraintSet{T}; title::String="",
    plt::Plots.Plot=plot(), save_figure::Bool=false) where {n1, T}

    n,m,N = size(solver)
    p  = get_model(solver).p
    pu = get_model(solver).pu
    px = get_model(solver).px

    X = TO.states(solver)
    U = TO.controls(solver)./2000.
    X_traj = [[[X[k][j] for k=1:N] for j in px[i]] for  i=1:p]
    U_traj = [[[U[k][j] for k=1:N-1] for j in pu[i]] for  i=1:p]
    colors = [:blue, :red, :yellow, :green, :green]

    # Boundary constraints
    for con in conSet.constraints
        if typeof(con.con) == BoundaryConstraint{T,2}
            i1 = findall(x->x==[con.con.x,con.con.y], px)[1]
            for j = 1:length(con.con)
                p1 = Array(con.con.p1)
                p2 = Array(con.con.p2)

                plot!([p1[1], p2[1]], [p1[2], p2[2]],
                    label="",
                    color=colors[i1],
                    marker=:circle,
                    linewidth=2.5)
            end
        end
    end

    # Circle constraints
    for con in conSet.constraints
        if typeof(con.con) in [CircleConstraint{T,length(con.con)}, ExpCircleConstraint{T,length(con.con)}]
            i1 = findall(x->x==[con.con.xi,con.con.yi], px)[1]
            for l in 1:length(con.con.x)
                plot!(circle_shape(con.con.x[l],
                                    con.con.y[l],
                                    con.con.radius[l]),
                    seriestype=[:shape,],
                    lw=2.5,
                    c=colors[i1],
                    linecolor=colors[i1],
                    label="",
                    fillalpha=0.3,
                    aspect_ratio=1)
            end
        end
    end

    # Plot arrows
    for i = 1:p
        for k = 1:N-1
            plot!(X_traj[i][1][k].+[0. ,U_traj[i][1][k]],
                X_traj[i][2][k].+[0. ,U_traj[i][2][k]],
                arrow = :arrow,
                linealpha = 0.5,
                linewidth = 4,
                linecolor = colors[i],
                label="")
        end
    end

    #Collision avoidance
    labels = ["Player $i" for i=1:p]
    for con in conSet.constraints
        if typeof(con.con) == CollisionConstraint{T,length(con.con)}
            i1 = findall(x->x==[con.con.x1,con.con.y1], px)[1]
            i2 = findall(x->x==[con.con.x2,con.con.y2], px)[1]
            for l in 1:length(con.con.x1)
                for k = 1:N
                    plot!(circle_shape(X_traj[i1][1][k],
                        X_traj[i1][2][k],
                        con.con.radius1[l]),
                        seriestype=[:shape,], lw=0.1, c=colors[i1],
                        linecolor=:black, label=labels[i1], fillalpha=0.2,
                        aspect_ratio=1)
                    plot!(circle_shape(X_traj[i2][1][k],
                        X_traj[i2][2][k],
                        con.con.radius2[l]),
                        seriestype=[:shape,], lw=0.1, c=colors[i2],
                        linecolor=:black, label=labels[i2], fillalpha=0.2,
                        aspect_ratio=1)
                        labels[i1] = ""
                        labels[i2] = ""
                end
            end
        end
    end

    # Ellipsoid Constraints
    actors_radii = get_actors_radii(solver)
    if any(actors_radii .!= actors_radii[1]*ones(p))
        @show "Inconsistent visualization due to different collision radii."
    end
    i_av = 1
    for con in solver.constraints.constraints
        if typeof(con.con) <: EllipsoidConstraint121
            inflation = actors_radii[i_av]*ones(length(px[i_av]))
            S = inv(inv(sqrt(con.con.S)) - Diagonal(inflation))
            for k in con.inds
                center = X[k][con.con.inds_human]
                plot_ellipse!(S, con.con.c, color=:green, M=200)
            end
        end
    end

    xlim = [Inf, -Inf]
    ylim = [Inf, -Inf]
    for i = 1:p
        xlim[1] = min(xlim[1], minimum(X_traj[i][1]))
        ylim[1] = min(ylim[1], minimum(X_traj[i][2]))
        xlim[2] = max(xlim[2], maximum(X_traj[i][1]))
        ylim[2] = max(ylim[2], maximum(X_traj[i][2]))
    end
    xlim += [-1.5, 1.5]
    ylim += [-0.5, 0.5]
    # Plot trajectory
    for i = 1:p
        plot!(X_traj[i][1], X_traj[i][2],
            title=title,
            color=:black,
            marker=:circle,
            linewidth=2.0,
            linestyle=:dot,
            xlim=xlim,
            ylim=ylim,
            label="")
        # Plot start and end
        plot!([x0[px[i][1]], xf[px[i][1]]],
            [x0[px[i][2]], xf[px[i][2]]],
            seriestype=:scatter,
            color=:orange,
            marker=:circle,
            linewidth=4.0,
            # label="StartEnd$i",
            label="")
    end
    if save_figure == true
        savefig("plots/car_trajectory_" * title * "_" * Dates.format(now(), "HH:MM:SS:sss") * ".jpg")
    end
    display(plt)
    return nothing
end

################################################################################
################################################################################
################################################################################

function get_actors_radii(solver::DirectGamesSolver{T}) where T
	model = get_model(solver)
	p = model.p
	px = model.px
	actors_radii = Vector{T}(undef, p)
	for con in solver.constraints.constraints
		if typeof(con.con) <: CollisionConstraint
			px1 = [con.con.x1, con.con.y1]
			px2 = [con.con.x2, con.con.y2]
			ind1 = findfirst(x -> x == px1, px)
			ind2 = findfirst(x -> x == px2, px)
			actors_radii[ind1] = con.con.radius1[1]
			actors_radii[ind2] = con.con.radius2[1]
		end
	end
	return actors_radii
end


function plot_ellipse!(S::AbstractMatrix{T}, c::AbstractArray{T}; M::Int=50, label::String="", color::Symbol=:blue) where T
	Θ = [(k-1)/(M-1)*2*pi for k=1:M]
	cir = [[cos(θ), sin(θ)] for θ in Θ]
	ell_points = [(cir[k] ./ norm(S*cir[k],2)) .+ c for k=1:M]
	ell_plot = [[e[1] for e in ell_points], [e[2] for e in ell_points]]
	plot!(ell_plot..., label=label, linewidth=0.5, color=color)
	# scatter!(ell_plot..., label=label, linewidth=0.5, color=color)
	return ell_points
end

################################################################################
################################################################################
################################################################################

function visualize_control(solver::TO.AbstractSolver{T}) where {T}
    visualize_control(TO.controls(solver), get_model(solver).pu)
end

function visualize_control(U::AbstractArray, pu::AbstractArray)
    plt = plot()
    p = length(pu)
    N = length(U)+1
    colors = [:blue, :red, :yellow, :green, :green]
    U_traj = [[[U[k][j] for k=1:N-1] for j in pu[i]] for  i=1:p]
    for i = 1:p
        first_plot = true
        for j = 1:length(pu[i])
            first_plot ? label = "player_$i" : label = ""
            first_plot = false
            plot!(U_traj[i][j],
                color=colors[i],
                marker=:circle,
                linewidth=1.0+i,
                linestyle=:dot,
                title="Control Inputs Trajectory",
                xlabel="Trajectory time steps",
                label=label)
        end
    end
    display(plt)
    return nothing
end

function visualize_state(solver::TO.AbstractSolver{T}) where {T}
    visualize_state(TO.states(solver))
end

function visualize_state(X::AbstractArray)
    plt = plot()
    N = length(X)
    n = length(X[1])
    colors = [:blue, :red, :yellow, :green, :green]
    first_plot = true
    x_traj = [[X[k][j] for k=1:N-1] for j=1:n]
    for j=1:n
        Plots.plot!(x_traj[j],
            color=colors[j%length(colors)+1],
            marker=:circle,
            linewidth=1.0+0.5*j,
            linestyle=:dot,
            title="State Trajectory",
            xlabel="Trajectory time steps",
            label="x_$j")
    end
    Plots.display(plt)
    return nothing
end

function visualize_collision_avoidance(solver::TO.AbstractSolver{T}) where T
    plt = plot()
    first_plot = true
    for (j,con) in enumerate(solver.constraints.constraints)
        if typeof(con.con) == CollisionConstraint{T,1}
            col = solver.constraints.constraints[j].vals
            λ_col = solver.constraints.constraints[j].λ
            col = [col[i][1] for i=1:length(col)]
            λ_col = [λ_col[i][1] for i=1:length(col)]

            first_plot ? labels = ["Collision constraint", "Collision multipliers"] : labels = ["",""]
            first_plot = false
            plot!(col,
                color=:orange,
                marker=:circle,
                linewidth=1.0,
                linestyle=:dot,
                title="Collision Avoidance Constraint",
                xlabel="Trajectory time steps",
                label=labels[1],
                legend=:bottomleft)
            plot!(λ_col,
                color=:blue,
                marker=:circle,
                linewidth=1.0,
                linestyle=:dot,
                label=labels[2])
        end
    end
    display(plt)
    return nothing
end

function visualize_circle_collision(solver::TO.AbstractSolver{T}) where T
    plt = plot()
    first_plot = true
    for con in solver.constraints.constraints
        if typeof(con.con) == CircleConstraint{T,length(con.con)}
            for j = 1:length(con.con)
                cir = con.vals
                λ_cir = con.λ
                cir = [cir[i][j] for i=1:length(cir)]
                λ_cir = [λ_cir[i][j] for i=1:length(cir)]

                first_plot ? labels = ["Collision constraint", "Collision multipliers"] : labels = ["",""]
                first_plot = false
                plot!(cir,
                    color=:orange,
                    marker=:circle,
                    linewidth=1.0,
                    linestyle=:dot,
                    title="Circle Collision Constraint",
                    xlabel="Trajectory time steps",
                    label=labels[1],
                    legend=:bottomleft,)
                plot!(λ_cir,
                    color=:blue,
                    marker=:circle,
                    linewidth=1.0,
                    linestyle=:dot,
                    label=labels[2])
            end
        end
    end
    display(plt)
    return nothing
end

function visualize_boundary_collision(solver::TO.AbstractSolver{T}) where T
    plt = plot()
    first_plot = true
    for con in solver.constraints.constraints
        if typeof(con.con) == BoundaryConstraint{T,2}
            for j = 1:length(con.con)
                bound = con.vals
                λ_bound = con.λ
                bound = [bound[i][j] for i=1:length(bound)]
                λ_bound = [λ_bound[i][j] for i=1:length(bound)]
                first_plot ? labels = ["Collision constraint", "Collision multipliers"] : labels = ["",""]
                first_plot = false
                plot!(bound,
                    color=:orange,
                    marker=:circle,
                    linewidth=1.0,
                    linestyle=:dot,
                    title="Boundary Collision Constraint",
                    xlabel="Trajectory time steps",
                    label=labels[1],
                    legend=:bottomleft,
                    )
                plot!(λ_bound,
                    color=:blue,
                    marker=:circle,
                    linewidth=1.0,
                    linestyle=:dot,
                    label=labels[2],
                    )
            end
        end
    end
    display(plt)
    return nothing
end

function visualize_dynamics(solver::TO.AbstractSolver{T}) where T
    plt = plot()
    n,m,pu,p = size(solver.model)
    dyn = solver.dyn_constraints.constraints[1].vals
    dyn = [log(10, mean(abs.(dyn[j]))) for j=1:length(dyn)]
    λ_dyn = [[mean(abs.(solver.ν[i][j])) for j=1:length(dyn)] for i = 1:p]
    plot!(dyn,
        color=:orange,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:left,
        title="Dynamics Constraint Violation and Multipliers",
        xlabel="Trajectory time steps",
        label="log( ||dynamics_viol.||_1)")
    plot!(λ_dyn,
        color=:blue,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:bottomleft,
        label="dynamics_multipliers")
    display(plt)
    return nothing
end

function visualize_optimality_merit(solver::TO.AbstractSolver)
    plt = plot()
    opt = zeros(solver.stats.iterations_total)
    count = 0
    for i = 1:solver.stats.iterations
        for j = 1:solver.stats.iterations_inner[i]
            count += 1
            opt[count] = log(10, solver.stats.optimality_merit[i][j])
        end
    end
    count = 0
    for i = 1:solver.stats.iterations
        inner_iter = solver.stats.iterations_inner[i]
        opt_f = filter(!isnan,opt)
        if length(opt_f) > 0
            maximum(opt_f)
            x = count
            y = minimum(opt_f)
            w = inner_iter
            h = maximum(opt_f) - minimum(opt_f) + 1.0
            plot!(rectangle_shape(x,y,w,h),
                label="Outer Loop $i",
                opacity=0.3)
        end
        count += inner_iter
    end
    plot!(opt,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        title="Optimality Constraint Satisfaction",
        ylabel="log(||G||_2)",
        xlabel="Solver Iterations",
        label="log(||G||_2)")
    display(plt)
    return nothing
end

function visualize_H_cond(solver::TO.AbstractSolver)
    plt = plot()
    cond = zeros(solver.stats.iterations_total)
    count = 0
    for i = 1:solver.stats.iterations
        for j = 1:solver.stats.iterations_inner[i]
            count += 1
            cond[count] = log(10, solver.stats.H_cond[i][j])
        end
    end
    count = 0
    for i = 1:solver.stats.iterations
        inner_iter = solver.stats.iterations_inner[i]
        cond_f = filter(!isnan,cond)
        if length(cond_f) > 0
            maximum(cond_f)
            x = count
            y = minimum(cond_f)
            w = inner_iter
            h = maximum(cond_f) - minimum(cond_f) + 1.0
            plot!(rectangle_shape(x,y,w,h),
                label="Outer Loop $i",
                opacity=0.3)
        end
        count += inner_iter
    end
    plot!(cond,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:topleft,
        title="Condition Number of Newton's Method",
        xlabel="Solver iterations",
        ylabel="log(cond(H))",
        label="log(cond(H))")
    display(plt)
    return nothing
end

function visualize_α(solver::TO.AbstractSolver)
    plt = plot()
    α = zeros(solver.stats.iterations_total)
    count = 0
    for i = 1:solver.stats.iterations
        for j = 1:solver.stats.iterations_inner[i]
            count += 1
            α[count] = log(2, solver.stats.α[i][j])
        end
    end
    count = 0
    for i = 1:solver.stats.iterations
        inner_iter = solver.stats.iterations_inner[i]
        α_f = filter(!isnan,α)
        if length(α_f) > 0
            maximum(α_f)
            x = count
            y = minimum(α_f)
            w = inner_iter
            h = maximum(α_f) - minimum(α_f) + 1.0
            plot!(rectangle_shape(x,y,w,h),
                label="Outer Loop $i",
                opacity=0.3)
        end
        count += inner_iter
    end
    plot!(α,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        label="log_2(alpha)",
        ylabel="log_2(alpha)",
        xlabel="Solver Iterations",
        title="Log_2 Newton Step Size")
    display(plt)
    return nothing
end

function visualize_α(solver::PenaltyiLQGamesSolver)
    plt = plot()
    α = log.(2, solver.stats.α)
    plot!(α,
        legend=false,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        ylabel="log_2(alpha)",
        xlabel="Solver Iterations",
        title="Log_2 Newton Step Size")
    display(plt)
    return nothing
end

function visualize_cmax(solver::TO.AbstractSolver)
    plt = plot()
    cmax = log.(10, solver.stats.cmax)
    plot!(cmax,
        legend=false,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        ylabel="log|Constraint_violation|_inf",
        xlabel="Solver Iterations",
        title="Log Maximum Constraint Violation")
    display(plt)
    return nothing
end

function circle_shape(h,k,r)
    θ = LinRange(0, 2*π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

function rectangle_shape(x,y,w,h)
    return Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
end
