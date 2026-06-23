"""
    plot_topology(flow; resolution=41, m_boundary=200, seed_count=20)

Visualize:
- vector field
- critical points
- boundary segments
- separatrix seeds
"""
function plot_topology(flow::VCFlowData.InterpolatedFlow;
    resolution::Int=41,
    zero_cell_tol=1e-10,
    duplicate_tol=1e-3,
    ignore_masked_cells=true,
    ignore_boundary_points=true,
    min_mask_boundary_distance_cells::Real=0,
    patch_boundary=false,
    critical_point_kwargs = (;)
)
    xmin, ymin, xmax, ymax = _spatial_bounds(flow)

    fig = Figure(size=(900, 700))
    ax = Axis(fig[1, 1],
        title="Vector Field",
        xlabel="x",
        ylabel="y",
        aspect=DataAspect()
    )

    # --- legend handles ---
    h_source  = scatter!(ax, [Point2f(0,0)]; color=:red, markersize=18, visible=false)
    h_sink    = scatter!(ax, [Point2f(0,0)]; color=:blue, markersize=18, visible=false)
    h_saddle  = scatter!(ax, [Point2f(0,0)]; color=:orange, marker=:diamond, markersize=20, visible=false)
    h_center  = scatter!(ax, [Point2f(0,0)]; color=:green, markersize=18, visible=false)
    h_sp_source  = scatter!(ax, [Point2f(0,0)]; color=:magenta, markersize=18, visible=false)
    h_sp_sink = scatter!(ax, [Point2f(0,0)]; color=:purple, markersize=18, visible=false)
    h_bsp     = scatter!(ax, [Point2f(0,0)]; color=:lightgray, markersize=14, visible=false)

    xs = collect(range(xmin, xmax; length=resolution))
    ys = collect(range(ymin, ymax; length=resolution))

    pts = Point2f[]
    vecs = Vec2f[]

    for x in xs, y in ys
        v = _flow_value(flow, SVector{2,Float64}(x, y))
        push!(pts, Point2f(x, y))
        push!(vecs, Vec2f(v[1], v[2]))
    end

    arrows2d!(ax, pts, vecs; lengthscale=0.08, alpha=0.5)

    cps = critical_points(flow;
        zero_cell_tol=zero_cell_tol,
        duplicate_tol=duplicate_tol,
        ignore_masked_cells=ignore_masked_cells,
        ignore_boundary_points=ignore_boundary_points,
        min_mask_boundary_distance_cells=min_mask_boundary_distance_cells
    )

    for cp in cps
        p = Point2f(cp.x[1], cp.x[2])

        if cp.kind isa Source
            scatter!(ax, [p]; markersize=18, color=:red, critical_point_kwargs...)
        elseif cp.kind isa Sink
            scatter!(ax, [p]; markersize=18, color=:blue, critical_point_kwargs...)
        elseif cp.kind isa Saddle
            scatter!(ax, [p]; markersize=20, color=:orange, marker=:diamond, critical_point_kwargs...)
        elseif cp.kind isa SpiralSource
            scatter!(ax, [p]; markersize=18, color=:magenta, critical_point_kwargs...)
        elseif cp.kind isa SpiralSink
            scatter!(ax, [p]; markersize=18, color=:purple, critical_point_kwargs...)
        else
            scatter!(ax, [p]; markersize=18, color=:green, critical_point_kwargs...)
        end
    end

    segs = boundary_segments(flow)

    for seg in segs
        mid = 0.5 * (seg.p0 + seg.p1)
        beh = boundary_behavior(flow, mid, seg.normal)

        color =
            beh == :inflow  ? :purple :
            beh == :outflow ? :darkgreen :
                              :gray

        lines!(ax, [seg.p0[1], seg.p1[1]], [seg.p0[2], seg.p1[2]]; color=color, linewidth=3)
    end

    bsps = boundary_switch_points(flow;
    patch=patch_boundary,
    include_mask_boundary=true,
    zero_cell_tol=zero_cell_tol
)

    for bsp in bsps
        scatter!(ax, [Point2f(bsp.x[1], bsp.x[2])]; markersize=14, color=:lightgray, critical_point_kwargs...)
    end

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)

    Legend(fig[1, 2],
        [h_source, h_sink, h_saddle, h_center, h_sp_source, h_sp_sink, h_bsp],
        ["Source", "Sink", "Saddle", "Center", "Spiral Source", "Spiral Sink", "Boundary Switch Point"]
    )

    return fig
end

"""
    plot_skeleton(flow; m_boundary=300)

Plot the topological skeleton:
- critical points
- boundary switch points
- separatrices between them
"""
function plot_skeleton(flow::VCFlowData.InterpolatedFlow;
    zero_cell_tol=1e-10,
    duplicate_tol=1e-3,
    ignore_masked_cells=true,
    ignore_boundary_points=true,
    min_mask_boundary_distance_cells::Real=0,
    patch_boundary=false,
    critical_point_kwargs = (;),
    separatrix_kwargs = (;),
)
    xmin, ymin, xmax, ymax = _spatial_bounds(flow)

    fig = Figure(size=(1000, 700))
    ax = Axis(fig[1, 1],
        title="Topological Skeleton",
        xlabel="x",
        ylabel="y",
        aspect=DataAspect()
    )

    h_source = scatter!(ax, [Point2f(0,0)]; color=:red, markersize=18, visible=false)
    h_sink   = scatter!(ax, [Point2f(0,0)]; color=:blue, markersize=18, visible=false)
    h_saddle = scatter!(ax, [Point2f(0,0)]; color=:yellow, marker=:circle, markersize=20, visible=false)
    h_center = scatter!(ax, [Point2f(0,0)]; color=:green, markersize=18, visible=false)
    h_sp_source = scatter!(ax, [Point2f(0,0)]; color=:magenta, markersize=18, visible=false)
    h_sp_sink = scatter!(ax, [Point2f(0,0)]; color=:purple, markersize=18, visible=false)
    h_bsp    = scatter!(ax, [Point2f(0,0)]; color=:lightgray, markersize=14, visible=false)


    cps = critical_points(flow;
        zero_cell_tol=zero_cell_tol,
        duplicate_tol=duplicate_tol,
        ignore_masked_cells=ignore_masked_cells,
        ignore_boundary_points=ignore_boundary_points,
        min_mask_boundary_distance_cells=min_mask_boundary_distance_cells
    )

    bsps = boundary_switch_points(flow;
    patch=patch_boundary,
    include_mask_boundary=true,
    zero_cell_tol=zero_cell_tol
)

    # boundary box 
    lines!(ax,
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin];
        color=:orange,
        linewidth=2
    )

    # separatrices from saddles 
    for cp in cps
        if cp.kind isa Saddle
            println("Tracing saddle at ", cp.x)

            for (x0, dir) in separatrix_seeds(flow, cp; ϵ=5e-3)
                pts = trace_separatrix(flow, x0;
                    dir=dir,
                    h=0.005,
                    stop_eps=5e-3,
                    minsteps_before_stop=10
                )

                println("  seed=", x0, " dir=", dir, " points=", length(pts))

                if length(pts) >= 2
                    lines!(ax, first.(pts), last.(pts); color=:black, linewidth=2)
                end
            end
        end
    end

    # separatrices from boundary switch points
    for bsp in bsps
        for (x0, dir) in separatrix_seeds(flow, bsp; ϵ=1e-3)
            pts = trace_separatrix(flow, x0;
                dir=dir,
                h=0.005,
                stop_eps=5e-3,
                minsteps_before_stop=10
            )

            println("BSP seed=", x0, " dir=", dir, " points=", length(pts))

            if length(pts) >= 2
                lines!(ax, first.(pts), last.(pts); color=:black, linewidth=2, separatrix_kwargs...)
            end
        end
    end

    println("Critical points:")
    for cp in cps
        println(cp.x, "  ", typeof(cp.kind))
    end

println("Boundary switch points: ", length(bsps))

    # critical points
    for cp in cps
        p = Point2f(cp.x[1], cp.x[2])

        if cp.kind isa Source
            scatter!(ax, [p]; markersize=18, color=:red, critical_point_kwargs...)
        elseif cp.kind isa Sink
            scatter!(ax, [p]; markersize=18, color=:blue, critical_point_kwargs...)
        elseif cp.kind isa Saddle
            scatter!(ax, [p]; markersize=20, color=:yellow, critical_point_kwargs...)
        elseif cp.kind isa SpiralSource
            scatter!(ax, [p]; markersize=18, color=:magenta, critical_point_kwargs...)
        elseif cp.kind isa SpiralSink
            scatter!(ax, [p]; markersize=18, color=:purple, critical_point_kwargs...)
        else
            scatter!(ax, [p]; markersize=18, color=:green, critical_point_kwargs...)
        end
    end

    # boundary switch points 
    for bsp in bsps
        scatter!(ax, [Point2f(bsp.x[1], bsp.x[2])];
            markersize=14,
            color=:lightgray,
            critical_point_kwargs...
        )
    end

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)

    Legend(fig[1, 2],
        [
            h_source,
            h_sink,
            h_saddle,
            h_center,
            h_sp_source,
            h_sp_sink,
            h_bsp
        ],
        [
            "Source",
            "Sink",
            "Saddle",
            "Center",
            "Spiral Source",
            "Spiral Sink",
            "Boundary Switch Point"
        ]
    )

    return fig
end

export plot_topology, plot_skeleton