using StaticArrays
using LinearAlgebra
using VCFlowData
using CairoMakie
using RK43


# Evaluate a stationary 2D VCFlowData.InterpolatedFlow at x
function _flow_value(flow::VCFlowData.InterpolatedFlow, x::SVector{2,T}) where {T}
    return flow.itp(x[1], x[2])
end

# Extract 2D spatial bounds
function _spatial_bounds(flow::VCFlowData.InterpolatedFlow)
    xmin, ymin = flow.lo
    xmax, ymax = flow.hi
    return xmin, ymin, xmax, ymax
end

function _inside(flow::VCFlowData.InterpolatedFlow, x::SVector{2})
    xmin, ymin, xmax, ymax = _spatial_bounds(flow)
    return xmin <= x[1] <= xmax && ymin <= x[2] <= ymax
end

function _safe_normalize(v::SVector{2,T}) where {T}
    n = norm(v)
    return n == 0 ? v : v / n
end

_linspace(a::T, b::T, n::Int) where {T} =
    n == 1 ? T[a] : collect(range(a, b; length=n))

function _classify_eigenvalues(λ; tol=1e-10)
    re = real.(λ)
    im = imag.(λ)

    # --- check for swirling ---
    is_complex = any(abs.(im) .> tol)

    if is_complex
        # complex conjugate pair -> swirling
        if all(re .> tol)
            return SpiralSource()
        elseif all(re .< -tol)
            return SpiralSink()
        else
            return Center()   # purely imaginary -> neutral swirl
        end
    else
        # purely real eigenvalues -> non-swirling
        if all(re .> tol)
            return Source()
        elseif all(re .< -tol)
            return Sink()
        elseif any(re .> tol) && any(re .< -tol)
            return Saddle()
        else
            return Center()
        end
    end
end

_tangent_from_normal(n::SVector{2,T}) where {T} = SVector{2,T}(-n[2], n[1])

function _signed_normal_component(flow::VCFlowData.InterpolatedFlow, x, normal)
    v = _flow_value(flow, x)
    return dot(v, normal)
end

"""
    jacobian(flow, t, x)

Jacobian matrix of the vector field at `(t, x)`.
For stationary flows, `t` is ignored.
"""
function jacobian(flow::VCFlowData.InterpolatedFlow, x::SVector{2,T}; h::T = sqrt(eps(T))) where {T}
    J = zeros(T, 2, 2)

    for i in 1:2
        dx = ntuple(j -> j == i ? h : zero(T), 2)
        δ = SVector{2,T}(dx)

        fp = _flow_value(flow, x + δ)
        fm = _flow_value(flow, x - δ)

        J[:, i] = (fp - fm) / (2h)
    end

    return SMatrix{2,2,T}(J)
end

"""
    critical_points(flow; nx=40, ny=40, tol=1e-8, maxiter=25)

Find critical points by:
1. scanning a coarse grid for candidate cells,
2. refining with Newton iterations,
3. classifying by eigenvalues of the Jacobian.
"""
function critical_points(flow::VCFlowData.InterpolatedFlow; tol=1e-8, maxiter::Int=25)
    itp = flow.itp
    xmin, ymin, xmax, ymax = _spatial_bounds(flow)

    xs = collect(range(xmin, xmax; length=length(axes(itp, 1))))
    ys = collect(range(ymin, ymax; length=length(axes(itp, 2))))

    candidates = SVector{2,Float64}[]

    # coarse scan: sign change in both components inside a cell
    for i in 1:(length(xs)-1), j in 1:(length(ys)-1)
        x0, x1 = xs[i], xs[i+1]
        y0, y1 = ys[j], ys[j+1]

        v00 = itp(x0, y0)
        v10 = itp(x1, y0)
        v01 = itp(x0, y1)
        v11 = itp(x1, y1)

        us = (v00[1], v10[1], v01[1], v11[1])
        vs = (v00[2], v10[2], v01[2], v11[2])

        if (minimum(us) <= 0 <= maximum(us)) &&
           (minimum(vs) <= 0 <= maximum(vs))

            xm = (x0 + x1) / 2
            ym = (y0 + y1) / 2

            push!(candidates, SVector(xm, ym))
        end
    end

    cps = CriticalPoint{Float64,2}[]

    for xstart in candidates
        x = xstart
        ok = true

        for _ in 1:maxiter
            f = _flow_value(flow, x)

            if norm(f) < tol
                break
            end

            J = jacobian(flow, x; h=sqrt(eps(Float64)))
            cond(J) > 1e12 && (ok = false; break)

            Δ = J \ f

            a = 1.0
            while !_inside(flow, x - a * Δ) || norm(_flow_value(flow, x - a * Δ)) > norm(f)
                a *= 0.5
                a < 1e-6 && (ok = false; break)
            end            
            
            x = x - a * Δ

            _inside(flow, x) || (ok = false; break)
        end

        ok || continue
        norm(_flow_value(flow, x)) < tol || continue

        duplicate = any(cp -> norm(cp.x - x) < 1e-5, cps)
        duplicate && continue

        J = jacobian(flow, x; h=sqrt(eps(Float64)))
        λ = eigvals(Matrix(J))
        kind = _classify_eigenvalues(λ)

        push!(cps, CriticalPoint(x, kind))
    end

    return cps
end

critical_type(cp::CriticalPoint) = cp.kind

"""
    boundary_behavior(flow, t, x, normal)

Classify the boundary behavior at x using the outward normal.
Returns :inflow, :outflow, or :tangent.
"""
function boundary_behavior(flow::VCFlowData.InterpolatedFlow, x, normal; tol=1e-10)
    v = _flow_value(flow, x)
    s = dot(v, normal)

    s < -tol && return :inflow
    s >  tol && return :outflow
    return :tangent
end

function boundary_behavior(flow::VCFlowData.InterpolatedFlow, t, x, normal; tol=1e-10)
    v = _flow_value(flow, x)
    s = dot(v, normal)

    s < -tol && return :inflow
    s >  tol && return :outflow
    return :tangent
end

"""
    boundary_segments(flow; m=200)

Split the rectangular boundary into segments of uniform boundary behavior.
"""
function boundary_segments(flow::VCFlowData.InterpolatedFlow)
    itp = flow.itp
    xsamples = length(axes(itp, 1))
    ysamples = length(axes(itp, 2))
    xmin, ymin, xmax, ymax = _spatial_bounds(flow)
    T = Float64

    edges = (
        (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmin, y), SVector{2,T}(-1, 0)), # left
        (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmax, y), SVector{2,T}( 1, 0)), # right
        (_linspace(xmin, xmax, xsamples), x -> SVector{2,T}(x, ymin), SVector{2,T}( 0,-1)), # bottom
        (_linspace(xmin, xmax, xsamples), x -> SVector{2,T}(x, ymax), SVector{2,T}( 0, 1)), # top
    )

    segs = BoundarySegment{Float64,2}[]

    for (params, mkpt, normal) in edges
        pts = [mkpt(p) for p in params]
        labels = [boundary_behavior(flow, p, normal) for p in pts]

        run_start = 1
        for i in 2:(length(labels)+1)
            if i == length(labels)+1 || labels[i] != labels[run_start]
                push!(segs, BoundarySegment(pts[run_start], pts[i-1], normal))
                run_start = i
            end
        end
    end

    return segs
end

"""
    boundary_switch_points(flow; m=400, tol=1e-10)

Find boundary switch points by detecting sign changes of v·n along each edge.
"""
function boundary_switch_points(flow::VCFlowData.InterpolatedFlow; tol=1e-10, patch::Bool=false)
    itp = flow.itp
    xsamples = length(axes(itp, 1))
    ysamples = length(axes(itp, 2))

    xmin, ymin, xmax, ymax = _spatial_bounds(flow)
    T = Float64

    pts = BoundarySwitchPoint{Float64,2}[]

    if !patch
        edges = (
            (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmin, y), SVector{2,T}(-1, 0)), # left
            (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmax, y), SVector{2,T}( 1, 0)), # right
            (_linspace(xmin, xmax, xsamples), x -> SVector{2,T}(x, ymin), SVector{2,T}( 0,-1)), # bottom
            (_linspace(xmin, xmax, xsamples), x -> SVector{2,T}(x, ymax), SVector{2,T}( 0, 1)), # top
        )

        for (params, mkpt, normal) in edges
            tangent = _tangent_from_normal(normal)

            xs = [mkpt(p) for p in params]
            svals = [_signed_normal_component(flow, x, normal) for x in xs]

            for i in 1:(length(xs)-1)
                x0 = xs[i]
                x1 = xs[i+1]

                s0 = svals[i]
                s1 = svals[i+1]

                if abs(s0) < tol
                    push!(pts, BoundarySwitchPoint(x0, normal, tangent))
                    continue
                end

                if abs(s1) < tol
                    push!(pts, BoundarySwitchPoint(x1, normal, tangent))
                    continue
                end

                sign(s0) == sign(s1) && continue

                α = clamp(s0 / (s0 - s1), 0.0, 1.0)
                xsw = (1 - α) * x0 + α * x1

                push!(pts, BoundarySwitchPoint(xsw, normal, tangent))
            end
        end
    else
        # Patch case:
        # The outer top/bottom layers have modified normal components.
        # Therefore use the first inner sample line and project the flow
        # onto the boundary tangent.
        xparams = _linspace(xmin, xmax, xsamples)
        yparams = _linspace(ymin, ymax, ysamples)

        y_bottom_inner = yparams[2]
        y_top_inner = yparams[end - 1]

        edges = (
            (
                xparams,
                x -> SVector{2,T}(x, y_bottom_inner), # offset sample line
                x -> SVector{2,T}(x, ymin),           # projected boundary point
                SVector{2,T}(0, -1)                   # bottom normal
            ),
            (
                xparams,
                x -> SVector{2,T}(x, y_top_inner),    # offset sample line
                x -> SVector{2,T}(x, ymax),           # projected boundary point
                SVector{2,T}(0, 1)                    # top normal
            ),
        )

        for (params, offset_mkpt, boundary_mkpt, normal) in edges
            tangent = _tangent_from_normal(normal)

            xs_offset = [offset_mkpt(p) for p in params]
            xs_boundary = [boundary_mkpt(p) for p in params]

            # 1D vector field along the offset curve
            svals = [dot(_flow_value(flow, x), tangent) for x in xs_offset]

            for i in 1:(length(xs_offset)-1)
                x0 = xs_boundary[i]
                x1 = xs_boundary[i+1]

                s0 = svals[i]
                s1 = svals[i+1]

                if abs(s0) < tol
                    push!(pts, BoundarySwitchPoint(x0, normal, tangent))
                    continue
                end

                if abs(s1) < tol
                    push!(pts, BoundarySwitchPoint(x1, normal, tangent))
                    continue
                end

                sign(s0) == sign(s1) && continue

                α = clamp(s0 / (s0 - s1), 0.0, 1.0)
                xsw = (1 - α) * x0 + α * x1

                push!(pts, BoundarySwitchPoint(xsw, normal, tangent))
            end
        end
    end

    out = BoundarySwitchPoint{Float64,2}[]

    for p in pts
        duplicate = any(q -> norm(q.x - p.x) < 1e-4, out)
        duplicate || push!(out, p)
    end

    return out
end

"""
    separatrix_seeds(flow, cp::CriticalPoint; ϵ=1e-6)

For saddles:
- unstable directions -> :forward
- stable directions   -> :backward
"""
function separatrix_seeds(flow::VCFlowData.InterpolatedFlow, cp::CriticalPoint; ϵ=1e-6)
    cp.kind isa Saddle || return Tuple{SVector{2,Float64},Symbol}[]

    x0 = SVector{2,Float64}(cp.x)
    J = jacobian(flow, x0)
    F = eigen(Matrix(J))

    seeds = Tuple{SVector{2,Float64},Symbol}[]

    for i in eachindex(F.values)
        λ = F.values[i]
        v = SVector{2,Float64}(F.vectors[:, i])
        v = _safe_normalize(v)
        v == SVector(0.0, 0.0) && continue

        if real(λ) > 0
            push!(seeds, (x0 + ϵ*v, :forward))
            push!(seeds, (x0 - ϵ*v, :forward))
        elseif real(λ) < 0
            push!(seeds, (x0 + ϵ*v, :backward))
            push!(seeds, (x0 - ϵ*v, :backward))
        end
    end

    return seeds
end

"""
    separatrix_seeds(flow, seg::BoundarySegment; k=20, ϵ=1e-6)

Sample inflow boundary segments and move seeds slightly into the domain.
"""
function separatrix_seeds(flow::VCFlowData.InterpolatedFlow, seg::BoundarySegment; k=20, ϵ=1e-6)
    mid = 0.5 * (seg.p0 + seg.p1)
    beh = boundary_behavior(flow, mid, seg.normal)

    beh == :inflow || return Tuple{SVector{2,Float64},Symbol}[]

    n̂ = _safe_normalize(SVector{2,Float64}(seg.normal))
    seeds = Tuple{SVector{2,Float64},Symbol}[]

    for s in range(0.0, 1.0; length=k)
        x = (1-s) * SVector{2,Float64}(seg.p0) + s * SVector{2,Float64}(seg.p1)
        x_in = x - ϵ * n̂
        push!(seeds, (x_in, :forward))
    end

    return seeds
end

"""
    separatrix_seeds(flow, bsp::BoundarySwitchPoint; ϵ=1e-3)

Create one seed slightly inside the domain near a boundary switch point.
From this seed, the separatrix should be traced both forward and backward.
"""
function separatrix_seeds(
    flow::VCFlowData.InterpolatedFlow,
    bsp::BoundarySwitchPoint;
    ϵ=1e-3
)
    x0 = SVector{2,Float64}(bsp.x)
    n̂ = _safe_normalize(SVector{2,Float64}(bsp.normal))

    # move only slightly into the domain
    xin = x0 - ϵ * n̂

    return [
        (xin, :forward),
        (xin, :backward),
    ]
end

integration_direction(::BoundarySegment) = :forward

"""
    trace_separatrix(flow, x0; dir=:forward, h=0.005, maxsteps=4000,
                     stop_eps=5e-3, minsteps_before_stop=10)

Trace a separatrix as sampled pathline using RK43.
"""
function trace_separatrix(flow::VCFlowData.InterpolatedFlow, x0::SVector{2,Float64};
    dir::Symbol=:forward,
    h::Float64=0.005,
    maxsteps::Int=4000,
    stop_eps::Float64=5e-3,
    minsteps_before_stop::Int=10
)
    cps = critical_points(flow)
    bsps = boundary_switch_points(flow)

    s = dir === :forward ? 1.0 : -1.0

    function dy(t, y)
        _inside(flow, y) || return RK43.OutOfDomain
        return s * _flow_value(flow, y)
    end

    opts = RK43.options(Float64;
        rtol=1e-4,
        atol=1e-7,
        hmax=h,
        maxsteps=maxsteps
    )

    solver = RK43.rk43solver(SVector{2,Float64}, opts)

    t0 = 0.0
    t1 = h * maxsteps

    state = RK43.initialize!(solver, dy, t0, t1, x0)
    state == RK43.Failed && return SVector{2,Float64}[x0]

    pts = SVector{2,Float64}[x0]

    for k in 1:maxsteps
        state = RK43.step!(solver, dy)

        if state == RK43.AcceptStep || state == RK43.Ok
            _, xnew, _ = RK43.tentativepos(solver)
            RK43.commit!(solver)

            _inside(flow, xnew) || break

            if k > minsteps_before_stop
                for cp in cps
                    if norm(cp.x - xnew) < stop_eps
                        push!(pts, SVector{2,Float64}(cp.x))
                        return pts
                    end
                end

                for bsp in bsps
                    if norm(bsp.x - xnew) < stop_eps
                        push!(pts, SVector{2,Float64}(bsp.x))
                        return pts
                    end
                end
            end

            push!(pts, xnew)

            RK43.isdone(solver) && break
        else
            break
        end
    end

    return pts
end

"""
    trace_saddle_separatrix(flow, cp, x0; dir=:forward, h=0.005,
                            maxsteps=4000, stop_eps=1e-2, minsteps_before_stop=10)

Trace a separatrix starting from a saddle seed.

Delegates the actual integration to `trace_separatrix`.
"""
function trace_saddle_separatrix(
    flow::VCFlowData.InterpolatedFlow,
    cp::CriticalPoint,
    x0::SVector{2,Float64};
    dir::Symbol=:forward,
    h::Float64=0.005,
    maxsteps::Int=4000,
    stop_eps::Float64=1e-2,
    minsteps_before_stop::Int=10
)
    return trace_separatrix(flow, x0;
        dir=dir,
        h=h,
        maxsteps=maxsteps,
        stop_eps=stop_eps,
        minsteps_before_stop=minsteps_before_stop
    )
end

function integration_direction(cp::CriticalPoint)
    cp.kind isa Source && return :backward
    cp.kind isa Sink   && return :forward
    cp.kind isa Saddle && return :both
    return :forward
end

