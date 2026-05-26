using TopologicalSkeleton
using Test
using StaticArrays

include(joinpath(@__DIR__, "..", "src", "abstracttopology.jl"))
include(joinpath(@__DIR__, "..", "src", "interpolatedtopology.jl"))
include(joinpath(@__DIR__, "..", "src", "flows", "simpleflow.jl"))
include(joinpath(@__DIR__, "..", "src", "decoder.jl"))

# Testflow 
flow = simpleflow()

@testset "Topology: Critical Points" begin

    cps = critical_points(flow)

    @test length(cps) == 3

    # expected positions
    expected = [
        SVector(-1.0, 0.0),
        SVector( 0.0, 0.0),
        SVector( 1.0, 0.0)
    ]

    # check positions (tolerant)
    for e in expected
        found = any(cp -> norm(cp.x - e) < 1e-2, cps)
        @test found
    end

    # check types
    for cp in cps
        if norm(cp.x - SVector(-1.0, 0.0)) < 1e-2
            @test cp.kind isa Saddle
        elseif norm(cp.x - SVector(0.0, 0.0)) < 1e-2
            @test cp.kind isa Sink
        elseif norm(cp.x - SVector(1.0, 0.0)) < 1e-2
            @test cp.kind isa Source
        end
    end
end

@testset "Topology: Boundary Segments" begin
    segs = boundary_segments(flow)

    @test length(segs) > 0

    for seg in segs
        mid = 0.5 * (seg.p0 + seg.p1)
        beh = boundary_behavior(flow, mid, seg.normal)
        @test beh in (:inflow, :outflow, :tangent)
    end
end

@testset "Topology: Separatrix Seeds" begin
    cps = critical_points(flow)

    for cp in cps
        if cp.kind isa Saddle
            seeds = separatrix_seeds(flow, cp)
            @test length(seeds) ≥ 2
        end
    end
end

@testset "loadflow formula" begin
    formula(x, y) = @SVector [x^3 - x, (x - 0.5) * y]
    xmin, xmax, ymin, ymax = -2.0, 2.0, -2.0, 2.0
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
    xmin, xmax, ymin, ymax = -2.0, 2.0, -2.0, 2.0
    flow = loadflow(matrix, xmin, xmax, ymin, ymax, false)
    @test flow isa VCFlowData.InterpolatedFlow
    @test flow.lo == SVector(xmin, ymin)
    @test flow.hi == SVector(xmax, ymax)
    # rand1, rand2 = rand(1:401), rand(1:401)
    # @test flow.itp(rand1, rand2) == matrix[rand1, rand2]
end

@testset "loadflow file" begin
    flow = loadflow("pipedcylinder2d.nc", false)
    @test flow isa VCFlowData.InterpolatedFlow
    @test flow.lo == SVector(-0.5, -0.5)
    @test flow.hi == SVector(5.5, 1.5)
end