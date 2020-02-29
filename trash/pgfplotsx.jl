#
# axis_1 = @pgf Axis(
#     {
#         "grid=both",
#         "minor y tick num=1",
#         "yminorgrids=true",
#         "tick align=outside",
#         # "ylabel near ticks",
#         # "xlabel near ticks",
#         "x label style={at={(axis description cs:0.5,-0.08)},anchor=north}",
#         "y label style={at={(axis description cs:-0.08,0.5)},rotate=0,anchor=south}",
#         "xlabel={Unemployment}",
#         "ylabel={inflation}",
#         xmajorgrids = false,
#         xmin = 0.00,   xmax = 3.00,
#         ymin = 0.00,   ymax = 30.00,
#     },
#     Plot(
#         {
#             "ybar interval",
#             "mark=none",
#             "fill=blue!25",
#         },
#         Table(fit(Histogram, range(0; stop = 1, length = 100).^3, closed = :left))
#         ),
# )
#
#
# axis_2 = @pgf Axis(
#     {
#         "grid=both",
#         "minor y tick num=1",
#         "yminorgrids=true",
#         "tick align=outside",
#         xmajorgrids = false,
#         xmin = 0.00,   xmax = 3.00,
#         ymin = 0.00,   ymax = 30.00,
#     },
#     Plot(
#         {
#             "ybar interval",
#             "mark=none",
#             "fill=red!25",
#         },
#         Table(fit(Histogram, range(0; stop = 1, length = 100).^3, closed = :left))
#         ),
# )
#
#
# gp = @pgf gp = GroupPlot(
#     {
#     group_style = { group_size = "2 by 2", "horizontal sep=1.4cm", "vertical sep=1.4cm", },
#     # height = "6cm", #width = "6cm"
#     });
#
#
#
#
# push!(gp, axis_1)
# push!(gp, axis_2)
# push!(gp, axis_2)
# push!(gp, axis_1)
# gp
#
# pgfsave("my_tikz.tikz", gp; include_preamble=false, dpi = 600)
