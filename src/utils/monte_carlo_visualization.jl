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
        color=:red,
        marker=:circle,
        linewidth=1.0,
        legend=:bottomright,
        title="Log Maximum Constraint Violation",
        xlabel="log( || Constraint Violation ||_inf)",
        ylabel="Occurences",
        label="")
    display(plt)
    return nothing
end

function visualize_H_cond(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    H_cond = log.(10, sampler.stats.H_cond)
    histogram!(H_cond,
        bins=40,
        color=:orange,
        marker=:circle,
        linewidth=1.0,
        legend=:bottomright,
        title="Log Condition Number of Newton's Method",
        xlabel="log(cond(H))",
        ylabel="Occurences",
        label="")
    display(plt)
    return nothing
end

function visualize_solve_time(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    solve_time = sampler.stats.solve_time
    histogram!(solve_time,
        bins=40,
        color=:green,
        marker=:circle,
        linewidth=1.0,
        legend=:bottomright,
        title="Solve Time",
        xlabel="Solve Time in s",
        ylabel="Occurences",
        label="")
    display(plt)
    return nothing
end

function visualize_iterations_total(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    iterations_total = sampler.stats.iterations_total
    histogram!(iterations_total,
        bins=40,
        color=:blue,
        marker=:circle,
        linewidth=1.0,
        legend=:bottomright,
        title="Total Number of Iterations to Converge",
        xlabel="Total Number of Iterations",
        ylabel="Occurences",
        label="")
    display(plt)
    return nothing
end

function visualize_optimality_merit(sampler::MonteCarloSampler{T}; save=false) where T
    plt = plot()
    optimality_merit = log.(10, sampler.stats.optimality_merit)
    histogram!(optimality_merit,
        bins=40,
        color=:yellow,
        marker=:circle,
        linewidth=1.0,
        legend=:bottomright,
        title="Optimality Constraint Satisfaction",
        xlabel="log(||G||_2)",
        ylabel="Occurences",
        label="")
    display(plt)
    return nothing
end
