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
    xmin_, ymin_, xmax_, ymax_ = _spatialbounds(flow)

    @test xmin_ ≈ xmin
    @test xmax_ ≈ xmax
    @test ymin_ ≈ ymin
    @test ymax_ ≈ ymax

    @test _inside(flow, SVector(0.0, 0.0))
    @test !_inside(flow, SVector(10.0, 0.0))

    @test _safenormalize(SVector(3.0, 4.0)) ≈ SVector(0.6, 0.8)
    @test _safenormalize(SVector(0.0, 0.0)) == SVector(0.0, 0.0)
    

    ls = _linspace(0.0, 1.0, 5)
    @test length(ls) == 5
    @test ls[1] ≈ 0.0
    @test ls[end] ≈ 1.0

    v = _flowvalue(flow, SVector(0.0, 0.0))
    @test v isa SVector{2,Float64}
end

@testset "Topology: Eigenvalue Classification" begin
    @test _classifyeigenvalues([1.0, 2.0]) isa Source
    @test _classifyeigenvalues([-1.0, -2.0]) isa Sink
    @test _classifyeigenvalues([-1.0, 2.0]) isa Saddle
    @test _classifyeigenvalues([0.0, 0.0]) isa Center

    @test _classifyeigenvalues([1.0 + im, 1.0 - im]) isa SpiralSource
    @test _classifyeigenvalues([-1.0 + im, -1.0 - im]) isa SpiralSink
end

@testset "Topology: Critical Points" begin
    cps = criticalpoints(flow)

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
            @test criticaltype(cp) isa Saddle
        elseif norm(cp.x - SVector(0.0, 0.0)) < 1e-2
            @test cp.kind isa Sink
            @test criticaltype(cp) isa Sink
        elseif norm(cp.x - SVector(1.0, 0.0)) < 1e-2
            @test cp.kind isa Source
            @test criticaltype(cp) isa Source
        end
    end
end

@testset "Topology: Adjacency Matrix for Polynomial Flow" begin
    adj_formula(x, y) = @SVector [x^3 - x, (x - 0.5) * y]
    adj_flow = loadflow(adj_formula, -2.0, 2.0, -2.0, 2.0, 401, 401, false)

    cps = criticalpoints(adj_flow;
        ignore_masked_cells=false,
        ignore_boundary_points=false,
        duplicate_tol=1e-10
    )

    bsps = boundaryswitchpoints(adj_flow;
        patch=false,
        include_mask_boundary=false
    )

    @test length(cps) == 3
    @test length(bsps) == 2

    saddle_idx = findfirst(cp -> cp.kind isa Saddle && norm(cp.x - SVector(-1.0, 0.0)) < 1e-2, cps)
    sink_idx   = findfirst(cp -> cp.kind isa Sink   && norm(cp.x - SVector( 0.0, 0.0)) < 1e-2, cps)
    source_idx = findfirst(cp -> cp.kind isa Source && norm(cp.x - SVector( 1.0, 0.0)) < 1e-2, cps)

    top_bsp_idx = findfirst(bsp -> norm(bsp.x - SVector(0.5,  2.0)) < 1e-2, bsps)
    bot_bsp_idx = findfirst(bsp -> norm(bsp.x - SVector(0.5, -2.0)) < 1e-2, bsps)

    @test saddle_idx !== nothing
    @test sink_idx !== nothing
    @test source_idx !== nothing
    @test top_bsp_idx !== nothing
    @test bot_bsp_idx !== nothing

    # Fixed node order for this test:
    # 1 = saddle at (-1, 0)
    # 2 = sink at (0, 0)
    # 3 = source at (1, 0)
    # 4 = top boundary switch point at (0.5, 2)
    # 5 = bottom boundary switch point at (0.5, -2)

    function build_adjacency_matrix(n::Integer, edges; directed::Bool=false)
        A = falses(n, n)

        for (i, j) in edges
            if i < 1 || i > n || j < 1 || j > n
                throw(ArgumentError("edge ($i, $j) is outside valid node range 1:$n"))
            end

            A[i, j] = true

            if !directed
                A[j, i] = true
            end
        end

        return A
    end

    edges = [
        (1, 2), # saddle -- sink
        (2, 3), # sink -- source
        (3, 4), # source -- top BSP
        (3, 5), # source -- bottom BSP
    ]

    A = build_adjacency_matrix(5, edges; directed=false)

    expected = [
        false  true   false  false  false;
        true   false  true   false  false;
        false  true   false  true   true;
        false  false  true   false  false;
        false  false  true   false  false;
    ]

    @test A == expected
    @test size(A) == (5, 5)
    @test eltype(A) == Bool
    @test issymmetric(A)

    @test_throws ArgumentError build_adjacency_matrix(5, [(1, 6)]; directed=false)
    @test_throws ArgumentError build_adjacency_matrix(5, [(0, 2)]; directed=false)
end

@testset "Topology: Bilinear Critical Point" begin
    # Piecewise bilinear vector field on one cell
    #
    # v00 = ( 3/4,  3/4)
    # v01 = (-1/4, -9/4)
    # v10 = (-9/4, -1/4)
    # v11 = ( 3/4,  3/4)

    V = Array{SVector{2,Float64}}(undef, 2, 2)

    V[1, 1] = SVector( 3/4,  3/4)   # v00
    V[1, 2] = SVector(-1/4, -9/4)   # v01
    V[2, 1] = SVector(-9/4, -1/4)   # v10
    V[2, 2] = SVector( 3/4,  3/4)   # v11

    flow = loadflow(V, 0.0, 1.0, 0.0, 1.0, false)

    cps = criticalpoints(flow)

    @test length(cps) == 2

    # Expected critical points
    @test any(cp -> norm(cp.x - SVector(0.25, 0.25)) < 1e-8, cps)
    @test any(cp -> norm(cp.x - SVector(0.75, 0.75)) < 1e-8, cps)

    for cp in cps
        if norm(cp.x - SVector(0.25, 0.25)) < 1e-8
            @test cp.kind isa Sink
        elseif norm(cp.x - SVector(0.75, 0.75)) < 1e-8
            @test cp.kind isa Saddle
        else
            @test false
        end
    end
end

@testset "Topology: Boundary Behavior" begin
    xmin_, ymin_, xmax_, ymax_ = _spatialbounds(flow)

    @test boundarybehavior(flow, SVector(xmin_, 0.0), SVector(-1.0, 0.0)) in (:inflow, :outflow, :tangent)
    @test boundarybehavior(flow, SVector(xmax_, 0.0), SVector( 1.0, 0.0)) in (:inflow, :outflow, :tangent)
    @test boundarybehavior(flow, SVector(0.0, ymin_), SVector(0.0, -1.0)) in (:inflow, :outflow, :tangent)
    @test boundarybehavior(flow, SVector(0.0, ymax_), SVector(0.0,  1.0)) in (:inflow, :outflow, :tangent)
end

@testset "Topology: Boundary Segments" begin
    segs = boundarysegments(flow)

    @test length(segs) > 0

    for seg in segs
        @test seg isa BoundarySegment{Float64,2}
        @test norm(seg.normal) ≈ 1.0

        mid = 0.5 * (seg.p0 + seg.p1)
        beh = boundarybehavior(flow, mid, seg.normal)
        @test beh in (:inflow, :outflow, :tangent)
    end
end

@testset "Topology: Boundary Switch Points" begin
    bsps = boundaryswitchpoints(flow)

    @test bsps isa Vector{BoundarySwitchPoint{Float64,2}}

    for bsp in bsps
        @test bsp isa BoundarySwitchPoint{Float64,2}
        @test norm(bsp.normal) ≈ 1.0
        @test norm(bsp.tangent) ≈ 1.0
        @test abs(dot(bsp.normal, bsp.tangent)) < 1e-10
    end
end

@testset "Topology: Separatrix Seeds from Critical Points" begin
    cps = criticalpoints(flow)

    for cp in cps
        seeds = separatrixseeds(flow, cp)

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
    segs = boundarysegments(flow)

    for seg in segs
        seeds = separatrixseeds(flow, seg)

        @test seeds isa Vector{Tuple{SVector{2,Float64},Symbol}}

        for (x, dir) in seeds
            @test x isa SVector{2,Float64}
            @test dir == :forward
        end
    end
end

@testset "Topology: Separatrix Seeds from Boundary Switch Points" begin
    bsps = boundaryswitchpoints(flow)

    for bsp in bsps
        seeds = separatrixseeds(flow, bsp)

        @test seeds isa Vector{Tuple{SVector{2,Float64},Symbol}}
        @test length(seeds) ≥ 1

        for (x, dir) in seeds
            @test x isa SVector{2,Float64}
            @test dir in (:forward, :backward)
        end
    end
end

@testset "Topology: Integration Direction" begin
    cps = criticalpoints(flow)

    for cp in cps
        dir = integrationdirection(cp)
        @test dir in (:forward, :backward, :both)
    end

    seg = first(boundarysegments(flow))
    @test integrationdirection(seg) == :forward
end

@testset "Topology: Trace Separatrix" begin
    cps = criticalpoints(flow)
    saddle = first(filter(cp -> cp.kind isa Saddle, cps))
    seeds = separatrixseeds(flow, saddle)

    @test length(seeds) ≥ 2

    x0, dir = first(seeds)
    pts = traceseparatrix(flow, x0; dir=dir, h=0.005, maxsteps=100)

    @test pts isa Vector{SVector{2,Float64}}
    @test length(pts) ≥ 1
    @test pts[1] == x0
end

@testset "Topology: Trace Saddle Separatrix" begin
    cps = criticalpoints(flow)
    saddle = first(filter(cp -> cp.kind isa Saddle, cps))
    seeds = separatrixseeds(flow, saddle)

    x0, dir = first(seeds)
    pts = tracesaddleseparatrix(flow, saddle, x0; dir=dir, h=0.005, maxsteps=100)

    @test pts isa Vector{SVector{2,Float64}}
    @test length(pts) ≥ 1
    @test pts[1] == x0
end

@testset "Topology: Plot Functions" begin
    fig1 = plottopology(flow)
    fig2 = plotskeleton(flow)

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

@testset "divergence" begin
    formula(x, y) = @SVector [x^3 - x, (x - 0.5) * y]
    flow = loadflow(formula, xmin, xmax, ymin, ymax, 401, 401, false)

    randx = rand() * (xmax - xmin) + xmin
    randy = rand() * (ymax - ymin) + ymin
    div_val = divergence(flow, SVector(randx, randy))

    # Compute expected divergence analytically
    expected_div = 3 * randx^2 - 1 + (randx - 0.5)
    @test abs(div_val - expected_div) < 0.02
end

@testset "poincarereturn" begin
    formula(x, y) = @SVector [-y,x]
    flow = loadflow(formula, -5.0, 5.0, -5.0, 5.0, 1001, 1001)
    p0 = SVector(1.0, 0.0)
    t = SVector(0.0, 1.0)
    start = 0.5
    firstreturn = poincarereturn(flow, p0, t, start)
    secondreturn = poincarereturn(flow, p0, t, firstreturn)
    @test abs(start - secondreturn) < 0.001
end