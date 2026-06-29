using StaticArrays

abstract type  AbstractTopology <: VCFlowData.AbstractFlow end

_impl(topo::AbstractTopology) = _impl(parent(topo))

Base.parent(topo::AbstractTopology) = error(
    "Base.parent not implemented for $(typeof(topo))"
)

abstract type CriticalPointType end
struct Source   <: CriticalPointType end
struct Sink     <: CriticalPointType end
struct Saddle   <: CriticalPointType end
struct Center   <: CriticalPointType end
struct SpiralSource <: CriticalPointType end
struct SpiralSink   <: CriticalPointType end

struct CriticalPoint{T,N} <: AbstractTopology
    x::SVector{N,T}
    kind::CriticalPointType
end

struct BoundarySegment{T,N} <: AbstractTopology
    p0::SVector{N,T}
    p1::SVector{N,T}
    normal::SVector{N,T}
end

struct BoundarySwitchPoint{T,N} <: AbstractTopology
    x::SVector{N,T}
    normal::SVector{N,T}
    tangent::SVector{N,T}
end

"""
    criticalpoints(flow)

Return all critical points of the flow.
"""
function criticalpoints end

"""
    criticaltype(cp)

Return type of critical point (:source, :sink, :saddle, :center, :spiral_source, :spiral_sink).
"""
function criticaltype end

"""
    separatrixseeds(flow, feature; ϵ=1e-6)

Return seed points for separatrix integration.
"""
function separatrixseeds end

"""
    boundarysegments(flow)

Return geometric boundary segments of the spatial domain.
"""
function boundarysegments end

"""
    boundaryswitchpoints(flow)

Return boundary switch points, i.e. points on the boundary where the
normal component changes sign.
"""
function boundaryswitchpoints end

"""
    boundarybehavior(flow, t, x, normal)

Classify boundary behavior (:inflow, :outflow, :tangent).
"""
function boundarybehavior end

"""
    integrationdirection(feature)

Return :forward or :backward integration direction.
"""
function integrationdirection end