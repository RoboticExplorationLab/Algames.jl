export
    visualize_cmax,
    visualize_H_cond,
    visualize_solve_time,
    visualize_iterations_total,
    visualize_optimality_merit


function visualize_cmax(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    cmax = log.(10, sampler.stats.cmax)
    histogram!(cmax,
        bins=40,
        # color=:blue,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:bottomleft,
        title="Cmax",
        label="log cmax")
    display(plt)
    return nothing
end

function visualize_H_cond(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    H_cond = log.(10, sampler.stats.H_cond)
    histogram!(H_cond,
        bins=40,
        # color=:blue,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:bottomleft,
        title="H_cond",
        label="log H_cond")
    display(plt)
    return nothing
end

function visualize_solve_time(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    solve_time = sampler.stats.solve_time
    histogram!(solve_time,
        bins=40,
        # color=:blue,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:bottomleft,
        title="solve_time",
        label="solve_time")
    display(plt)
    return nothing
end

function visualize_iterations_total(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    iterations_total = sampler.stats.iterations_total
    histogram!(iterations_total,
        bins=40,
        # color=:blue,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:bottomleft,
        title="iterations_total",
        label="iterations_total")

    display(plt)
    return nothing
end

function visualize_optimality_merit(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    optimality_merit = log.(10, sampler.stats.optimality_merit)
    histogram!(optimality_merit,
        bins=40,
        # color=:blue,
        marker=:circle,
        linewidth=1.0,
        linestyle=:dot,
        legend=:bottomleft,
        title="optimality_merit",
        label="log optimality_merit")
    display(plt)
    return nothing
end
