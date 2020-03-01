export
    plot_velocity,
    latex_plot_velocity

function plot_velocity(solver_mpc, solver)
    plt = plot()
    idx = 4
    iter = solver_mpc.stats.iterations
    vel_real = [TO.states(solver_mpc)[k][idx] for k=1:iter]
    vel_pred = [TO.states(solver)[k][idx] for k=1:N]
    max_vel = max(maximum(vel_pred), maximum(vel_real))
    vel_real ./= max_vel
    vel_pred ./= max_vel


    time_real = cumsum(solver_mpc.stats.solve_time[1:iter]).- solver_mpc.stats.solve_time[1]
    time_pred = [(k-1)*solver.tf/(N-1) for k=1:N]
    max_time = max(maximum(time_pred), maximum(time_real))
    time_real ./= max_time
    time_pred ./= max_time
    plot!(time_pred, vel_pred,
        # color=:blue,
        # marker=:circle,
        linewidth=3.0,
        linestyle=:dot,
        title="Velocity",
        label="vel_pred")
    plot!(time_real, vel_real,
        # color=:blue,
        # marker=:circle,
        linewidth=3.0,
        linestyle=:dot,
        title="Velocity",
        label="vel_real")
    display(plt)
    return nothing
end

function latex_plot_velocity(solver_mpc, solver)
    idx = 4
    iter = solver_mpc.stats.iterations
    vel_real = [TO.states(solver_mpc)[k][idx] for k=1:iter]
    vel_pred = [TO.states(solver)[k][idx] for k=1:N]
    max_vel = max(maximum(vel_pred), maximum(vel_real))
    vel_real ./= max_vel
    vel_pred ./= max_vel


    time_real = cumsum(solver_mpc.stats.solve_time[1:iter]).- solver_mpc.stats.solve_time[1]
    time_pred = [(k-1)*solver.tf/(N-1) for k=1:N]
    max_time = max(maximum(time_pred), maximum(time_real))
    time_real ./= max_time
    time_pred ./= max_time

    x = range(-1; stop = 1, length = 51) # so that it contains 1/0
    axis = @pgf Axis(
        {
            "legend pos=south east",
            ymajorgrids,
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
            "x label style={at={(axis description cs:0.5,-0.20)},anchor=north}",
            "y label style={at={(axis description cs:-0.10,0.5)},rotate=0,anchor=south}",
            "xlabel={Scaled Time}",
            "ylabel={Scaled Velocity}",
            xmajorgrids = false,
            xmin = 0.00,   xmax = 1.00,
            ymin = 0.00, #  ymax = 1.00,

        },
        Plot(
            {
            "thick",
            "orange",
            no_marks,
            },
            Coordinates(time_pred, vel_pred)
        ),
        LegendEntry("Initial Plan"),
        Plot(
            {
            "thick",
            "blue",
            no_marks,
            },
            Coordinates(time_real, vel_real)
        ),
        LegendEntry("MPC"),
    )
    pgfsave("plots/tikz/velocity_"*string(solver_mpc.stats.iterations)*".tikz",
        axis; include_preamble=false, dpi = 600)
    return axis
end

# latex_plot_velocity(mpc_solver, solver_directgames_clean)
# plot_velocity(mpc_solver, solver_directgames_clean)
