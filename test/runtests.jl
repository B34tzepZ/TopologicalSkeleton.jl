using TopologicalSkeleton
using Test
using StaticArrays
using LinearAlgebra
using Interpolations
using NCDatasets
using VCDataSets
using VCFlowData
using CairoMakie

include(joinpath(@__DIR__, "..", "src", "abstracttopology.jl"))
include(joinpath(@__DIR__, "..", "src", "interpolatedtopology.jl"))
include(joinpath(@__DIR__, "..", "src", "decoder.jl"))
include(joinpath(@__DIR__, "..", "src", "plot.jl"))

formula(x, y) = @SVector [x^3 - x, (x - 0.5) * y]
xmin, xmax, ymin, ymax = -2.0, 2.0, -2.0, 2.0
flow = loadflow(formula, xmin, xmax, ymin, ymax, 401, 401, false)

@testset "Topology helpers" begin
    xmin_, ymin_, xmax_, ymax_ = _spatial_bounds(flow)

    @test xmin_ ≈ xmin
    @test xmax_ ≈ xmax
    @test ymin_ ≈ ymin
    @test ymax_ ≈ ymax

    @test _inside(flow, SVector(0.0, 0.0))
    @test !_inside(flow, SVector(10.0, 0.0))

    @test _safe_normalize(SVector(3.0, 4.0)) ≈ SVector(0.6, 0.8)
    @test _safe_normalize(SVector(0.0, 0.0)) == SVector(0.0, 0.0)
    

    ls = _linspace(0.0, 1.0, 5)
    @test length(ls) == 5
    @test ls[1] ≈ 0.0
    @test ls[end] ≈ 1.0

    v = _flow_value(flow, SVector(0.0, 0.0))
    @test v isa SVector{2,Float64}
end

@testset "Topology: Eigenvalue Classification" begin
    @test _classify_eigenvalues([1.0, 2.0]) isa Source
    @test _classify_eigenvalues([-1.0, -2.0]) isa Sink
    @test _classify_eigenvalues([-1.0, 2.0]) isa Saddle
    @test _classify_eigenvalues([0.0, 0.0]) isa Center

    @test _classify_eigenvalues([1.0 + im, 1.0 - im]) isa SpiralSource
    @test _classify_eigenvalues([-1.0 + im, -1.0 - im]) isa SpiralSink
end

@testset "Topology: Critical Points" begin
    cps = critical_points(flow)

    @test length(cps) == 3

    expected = [
        SVector(-1.0, 0.0),
        SVector( 0.0, 0.0),
        SVector( 1.0, 0.0)
    ]

    for e in expected
        found = any(cp -> norm(cp.x - e) < 1e-2, cps)
        @test found
    end

    for cp in cps
        if norm(cp.x - SVector(-1.0, 0.0)) < 1e-2
            @test cp.kind isa Saddle
            @test critical_type(cp) isa Saddle
        elseif norm(cp.x - SVector(0.0, 0.0)) < 1e-2
            @test cp.kind isa Sink
            @test critical_type(cp) isa Sink
        elseif norm(cp.x - SVector(1.0, 0.0)) < 1e-2
            @test cp.kind isa Source
            @test critical_type(cp) isa Source
        end
    end
end

@testset "Topology: Boundary Behavior" begin
    xmin_, ymin_, xmax_, ymax_ = _spatial_bounds(flow)

    @test boundary_behavior(flow, SVector(xmin_, 0.0), SVector(-1.0, 0.0)) in (:inflow, :outflow, :tangent)
    @test boundary_behavior(flow, SVector(xmax_, 0.0), SVector( 1.0, 0.0)) in (:inflow, :outflow, :tangent)
    @test boundary_behavior(flow, SVector(0.0, ymin_), SVector(0.0, -1.0)) in (:inflow, :outflow, :tangent)
    @test boundary_behavior(flow, SVector(0.0, ymax_), SVector(0.0,  1.0)) in (:inflow, :outflow, :tangent)
end

@testset "Topology: Boundary Segments" begin
    segs = boundary_segments(flow)

    @test length(segs) > 0

    for seg in segs
        @test seg isa BoundarySegment{Float64,2}
        @test norm(seg.normal) ≈ 1.0

        mid = 0.5 * (seg.p0 + seg.p1)
        beh = boundary_behavior(flow, mid, seg.normal)
        @test beh in (:inflow, :outflow, :tangent)
    end
end

@testset "Topology: Boundary Switch Points" begin
    bsps = boundary_switch_points(flow)

    @test bsps isa Vector{BoundarySwitchPoint{Float64,2}}

    for bsp in bsps
        @test bsp isa BoundarySwitchPoint{Float64,2}
        @test norm(bsp.normal) ≈ 1.0
        @test norm(bsp.tangent) ≈ 1.0
        @test abs(dot(bsp.normal, bsp.tangent)) < 1e-10
    end
end

@testset "Topology: Separatrix Seeds from Critical Points" begin
    cps = critical_points(flow)

    for cp in cps
        seeds = separatrix_seeds(flow, cp)

        @test seeds isa Vector{Tuple{SVector{2,Float64},Symbol}}

        if cp.kind isa Saddle
            @test length(seeds) ≥ 2
        else
            @test length(seeds) == 0
        end

        for (x, dir) in seeds
            @test x isa SVector{2,Float64}
            @test dir in (:forward, :backward)
        end
    end
end

@testset "Topology: Separatrix Seeds from Boundary Segments" begin
    segs = boundary_segments(flow)

    for seg in segs
        seeds = separatrix_seeds(flow, seg)

        @test seeds isa Vector{Tuple{SVector{2,Float64},Symbol}}

        for (x, dir) in seeds
            @test x isa SVector{2,Float64}
            @test dir == :forward
        end
    end
end

@testset "Topology: Separatrix Seeds from Boundary Switch Points" begin
    bsps = boundary_switch_points(flow)

    for bsp in bsps
        seeds = separatrix_seeds(flow, bsp)

        @test seeds isa Vector{Tuple{SVector{2,Float64},Symbol}}
        @test length(seeds) ≥ 1

        for (x, dir) in seeds
            @test x isa SVector{2,Float64}
            @test dir in (:forward, :backward)
        end
    end
end

@testset "Topology: Integration Direction" begin
    cps = critical_points(flow)

    for cp in cps
        dir = integration_direction(cp)
        @test dir in (:forward, :backward, :both)
    end

    seg = first(boundary_segments(flow))
    @test integration_direction(seg) == :forward
end

@testset "Topology: Trace Separatrix" begin
    cps = critical_points(flow)
    saddle = first(filter(cp -> cp.kind isa Saddle, cps))
    seeds = separatrix_seeds(flow, saddle)

    @test length(seeds) ≥ 2

    x0, dir = first(seeds)
    pts = trace_separatrix(flow, x0; dir=dir, h=0.005, maxsteps=100)

    @test pts isa Vector{SVector{2,Float64}}
    @test length(pts) ≥ 1
    @test pts[1] == x0
end

@testset "Topology: Trace Saddle Separatrix" begin
    cps = critical_points(flow)
    saddle = first(filter(cp -> cp.kind isa Saddle, cps))
    seeds = separatrix_seeds(flow, saddle)

    x0, dir = first(seeds)
    pts = trace_saddle_separatrix(flow, saddle, x0; dir=dir, h=0.005, maxsteps=100)

    @test pts isa Vector{SVector{2,Float64}}
    @test length(pts) ≥ 1
    @test pts[1] == x0
end

@testset "Topology: Plot Functions" begin
    fig1 = plot_topology(flow)
    fig2 = plot_skeleton(flow)

    @test fig1 isa Figure
    @test fig2 isa Figure
end

@testset "loadflow formula" begin
    formula(x, y) = @SVector [x^3 - x, (x - 0.5) * y]
    flow = loadflow(formula, xmin, xmax, ymin, ymax, 401, 401, false)

    @test flow isa VCFlowData.InterpolatedFlow
    @test flow.lo == SVector(xmin, ymin)
    @test flow.hi == SVector(xmax, ymax)

    randx = rand() * (xmax - xmin) + xmin
    randy = rand() * (ymax - ymin) + ymin
    @test all(abs.(flow.itp(randx, randy) - formula(randx, randy)) .< 0.001)
end

@testset "loadflow matrix" begin
    matrix = [@SVector rand(Float64, 2) for i in 1:401, j in 1:401]
    flow = loadflow(matrix, xmin, xmax, ymin, ymax, false)

    @test flow isa VCFlowData.InterpolatedFlow
    @test flow.lo == SVector(xmin, ymin)
    @test flow.hi == SVector(xmax, ymax)
end

@testset "loadflow file" begin
    flow = loadflow("pipedcylinder2d.nc", false)

    @test flow isa VCFlowData.InterpolatedFlow
    @test flow.lo == SVector(-0.5, -0.5)
    @test flow.hi == SVector(5.5, 1.5)
end