using StaticArrays

abstract type FlowFeature end

abstract type CriticalPointType end
struct Source   <: CriticalPointType end
struct Sink     <: CriticalPointType end
struct Saddle   <: CriticalPointType end
struct Center   <: CriticalPointType end
struct SpiralSource <: CriticalPointType end
struct SpiralSink   <: CriticalPointType end

struct CriticalPoint{T,N} <: FlowFeature
    x::SVector{N,T}
    kind::CriticalPointType
end

struct BoundarySegment{T,N} <: FlowFeature
    p0::SVector{N,T}
    p1::SVector{N,T}
    normal::SVector{N,T}
end

struct BoundarySwitchPoint{T,N} <: FlowFeature
    x::SVector{N,T}
    normal::SVector{N,T}
    tangent::SVector{N,T}
end

"""
    jacobian(flow, t, x)

Jacobian matrix of the vector field at (t, x).
"""
function jacobian end

"""
    critical_points(flow)

Return all critical points of the flow.
"""
function critical_points end

"""
    critical_type(cp)

Return type of critical point (:source, :sink, :saddle, :center, :spiral_source, :spiral_sink).
"""
function critical_type end

"""
    separatrix_seeds(flow, feature; ϵ=1e-6)

Return seed points for separatrix integration.
"""
function separatrix_seeds end

"""
    boundary_segments(flow)

Return geometric boundary segments of the spatial domain.
"""
function boundary_segments end

"""
    boundary_switch_points(flow)

Return boundary switch points, i.e. points on the boundary where the
normal component changes sign.
"""
function boundary_switch_points end

"""
    boundary_behavior(flow, t, x, normal)

Classify boundary behavior (:inflow, :outflow, :tangent).
"""
function boundary_behavior end

"""
    integration_direction(feature)

Return :forward or :backward integration direction.
"""
function integration_direction end