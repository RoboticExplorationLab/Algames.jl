export
	visualize_latex_cmax,
	visualize_latex_H_cond,
	visualize_latex_solve_time,
	visualize_latex_iterations_total,
	visualize_latex_optimality_merit,
	visualize_latex_sampler

function visualize_latex_cmax(sampler::MonteCarloSampler{T}; save=false) where T
    cmax = log.(10, sampler.stats.cmax)
    axis = @pgf Axis(
        {
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
			# "ymode=log",
            # "ylabel near ticks",
            # "xlabel near ticks",
            "x label style={at={(axis description cs:0.5,-0.18)},anchor=north}",
            "y label style={at={(axis description cs:-0.18,0.5)},rotate=0,anchor=south}",
            "xlabel={Log Constraint Violation}",
            # "ylabel={Occurrences}",
            xmajorgrids = false,
            xmin = -6.00,   xmax = 0.00,
            ymin =   0.00, #  ymax = 30.00,
        },
        Plot(
            {
                "ybar interval",
                "mark=none",
                "fill=blue!25",
            },
            Table(fit(Histogram, cmax, -6.00:0.30:0.00)) # , closed=:left))
            ),
    )
	if save
		@show "we save"
		pgfsave("plots/tikz/cmax_"*string(sampler.opts.iterations)*".tikz",
			axis; include_preamble=false, dpi = 600)
	end
    return axis
end

function visualize_latex_H_cond(sampler::MonteCarloSampler{T}; save=false) where T
    H_cond = log.(10, sampler.stats.H_cond)
    axis = @pgf Axis(
        {
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
			# "ymode=log",
            # "ylabel near ticks",
            # "xlabel near ticks",
            "x label style={at={(axis description cs:0.5,-0.18)},anchor=north}",
            "y label style={at={(axis description cs:-0.18,0.5)},rotate=0,anchor=south}",
            "xlabel={Log Condition Number}",
            # "ylabel={Occurrences}",
            xmajorgrids = false,
            xmin = 0.00,   xmax = 12.00,
            ymin = 0.00, #  ymax = 30.00,
        },
        Plot(
            {
                "ybar interval",
                "mark=none",
                "fill=red!25",
            },
            Table(fit(Histogram, H_cond, 0.00:0.60:12.00)) # , closed=:left))
            ),
    )
	if save
		@show "we save"
		pgfsave("plots/tikz/H_cond_"*string(sampler.opts.iterations)*".tikz",
			axis; include_preamble=false, dpi = 600)
	end
    return axis
end

function visualize_latex_solve_time(sampler::MonteCarloSampler{T}; save=false) where T
    solve_time = sampler.stats.solve_time
    axis = @pgf Axis(
        {
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
			# "ymode=log",
            # "ylabel near ticks",
            # "xlabel near ticks",
            "x label style={at={(axis description cs:0.5,-0.18)},anchor=north}",
            "y label style={at={(axis description cs:-0.18,0.5)},rotate=0,anchor=south}",
            "xlabel={Solve Time in s}",
            # "ylabel={Occurrences}",
            xmajorgrids = false,
            xmin = 0.00,   xmax = 6.00,
            ymin = 0.00, #  ymax = 30.00,
        },
        Plot(
            {
                "ybar interval",
                "mark=none",
                "fill=red!25",
            },
            Table(fit(Histogram, solve_time, 0.00:0.3:6.00)) # , closed=:left))
            ),
    )
	if save
		@show "we save"
		pgfsave("plots/tikz/solve_time_"*string(sampler.opts.iterations)*".tikz",
			axis; include_preamble=false, dpi = 600)
	end
    return axis
end

function visualize_latex_iterations_total(sampler::MonteCarloSampler{T}; save=false) where T
    iterations_total = sampler.stats.iterations_total
    axis = @pgf Axis(
        {
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
			# "ymode=log",
            # "ylabel near ticks",
            # "xlabel near ticks",
            "x label style={at={(axis description cs:0.5,-0.18)},anchor=north}",
            "y label style={at={(axis description cs:-0.18,0.5)},rotate=0,anchor=south}",
            "xlabel={Solver Iterations}",
            # "ylabel={Occurrences}",
            xmajorgrids = false,
            xmin = 0.00,   xmax = 160.00,
            ymin = 0.00, #  ymax = 30.00,
        },
        Plot(
            {
                "ybar interval",
                "mark=none",
                "fill=orange!50",
            },
            Table(fit(Histogram, iterations_total, 0.00:8.00:160.00)) # , closed=:left))
            ),
    )
	if save
		@show "we save"
		pgfsave("plots/tikz/iterations_"*string(sampler.opts.iterations)*".tikz",
			axis; include_preamble=false, dpi = 600)
	end
    return axis
end

function visualize_latex_optimality_merit(sampler::MonteCarloSampler{T}; save=false) where T
    optimality_merit = log.(10, sampler.stats.optimality_merit)
    axis = @pgf Axis(
        {
            "grid=both",
            "minor y tick num=1",
            "yminorgrids=true",
            "tick align=outside",
			# "ymode=log",
            # "ylabel near ticks",
            # "xlabel near ticks",
            "x label style={at={(axis description cs:0.5,-0.18)},anchor=north}",
            "y label style={at={(axis description cs:-0.18,0.5)},rotate=0,anchor=south}",
            "xlabel={Log Optimality Criterion}",
            # "ylabel={Occurrences}",
            xmajorgrids = false,
            xmin = -4.00,   xmax = 4.00,
            ymin = 0.00, #  ymax = 30.00,
        },
        Plot(
            {
                "ybar interval",
                "mark=none",
                "fill=green!50",
            },
            Table(fit(Histogram, optimality_merit, -4.00:0.40:4.00)) # , closed=:left))
            ),
    )
	if save
		@show "we save"
		pgfsave("plots/tikz/optimality_merit_"*string(sampler.opts.iterations)*".tikz",
			axis; include_preamble=false, dpi = 600)
	end
    return axis
end

function visualize_latex_sampler(sampler::MonteCarloSampler{T}; save=false) where T
	# axis_H_cond = visualize_latex_H_cond(sampler)
	axis_solve_time = visualize_latex_solve_time(sampler)
	axis_cmax = visualize_latex_cmax(sampler)
	axis_iterations_total = visualize_latex_iterations_total(sampler)
	axis_optimality_merit = visualize_latex_optimality_merit(sampler)

    if string(typeof(sampler.solver).name) == "DirectGamesSolver"
	    gp = @pgf gp = GroupPlot(
	        {
	        group_style = { group_size = "2 by 2", "horizontal sep=1.4cm", "vertical sep=1.4cm", },
	        # height = "6cm", #width = "6cm"
	        });
			push!(gp, axis_solve_time)
			push!(gp, axis_iterations_total)
			push!(gp, axis_cmax)
			push!(gp, axis_optimality_merit)
			if save
				pgfsave("plots/tikz/sampler_direct_"*string(sampler.opts.iterations)*".tikz",
					gp; include_preamble=false, dpi = 600)
			end
	elseif string(typeof(sampler.solver).name) == "PenaltyiLQGamesSolver"
	    gp = @pgf gp = GroupPlot(
	        {
	        group_style = { group_size = "2 by 2", "horizontal sep=1.4cm", "vertical sep=1.4cm", },
	        # height = "6cm", #width = "6cm"
	        });
		empty_axis = @pgf Axis()
	    push!(gp, axis_solve_time)
	    push!(gp, axis_iterations_total)
		push!(gp, axis_cmax)

		if save
			@show "we save"
	        pgfsave("plots/tikz/sampler_iLQGames_"*string(sampler.opts.iterations)*".tikz",
	            gp; include_preamble=false, dpi = 600)
	    end
	end
    return gp
end
