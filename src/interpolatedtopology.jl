using StaticArrays
using LinearAlgebra
using VCFlowData
using CairoMakie
using RK43
using ForwardDiff

# Evaluate a stationary 2D VCFlowData.InterpolatedFlow at x
function _flowvalue(flow::VCFlowData.InterpolatedFlow, x::SVector{2,T}) where {T}
    return flow.itp(x[1], x[2])
end

# Extract 2D spatial bounds
function _spatialbounds(flow::VCFlowData.InterpolatedFlow)
    xmin, ymin = flow.lo
    xmax, ymax = flow.hi
    return xmin, ymin, xmax, ymax
end

function _inside(flow::VCFlowData.InterpolatedFlow, x::SVector{2})
    xmin, ymin, xmax, ymax = _spatialbounds(flow)
    return xmin <= x[1] <= xmax && ymin <= x[2] <= ymax
end

function _safenormalize(v::SVector{2,T}) where {T}
    n = norm(v)
    return n == 0 ? v : v / n
end

_linspace(a::T, b::T, n::Int) where {T} =
    n == 1 ? T[a] : collect(range(a, b; length=n))

function _classifyeigenvalues(λ; tol=1e-10)
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

function divergence(flow::VCFlowData.InterpolatedFlow, x::SVector{2,T}) where {T}
    J = ForwardDiff.jacobian(x -> _flowvalue(flow, x), x)
    return tr(J)
end

function _pointsegmentdistance(p::SVector{2,Float64}, a::SVector{2,Float64}, b::SVector{2,Float64})
    ab = b - a
    denom = dot(ab, ab)

    denom == 0 && return norm(p - a)

    t = clamp(dot(p - a, ab) / denom, 0.0, 1.0)
    q = a + t * ab

    return norm(p - q)
end

"""
    criticalpoints(flow; tol=1e-8)

Find critical points of a 2D bilinearly interpolated vector field by:

1. constructing the bilinear interpolation inside every grid cell,
2. solving the resulting bilinear system analytically,
3. locating all interior zeros of the vector field,
4. classifying them using the eigenvalues of the local Jacobian.
"""
function criticalpoints(flow::VCFlowData.InterpolatedFlow;
    tol=1e-8,
    zero_cell_tol=1e-12,
    duplicate_tol=1e-4,
    ignore_masked_cells::Bool=true,
    ignore_boundary_points::Bool=true,
    min_mask_boundary_distance_cells::Real=0,
)
    itp = flow.itp
    xmin, ymin, xmax, ymax = _spatialbounds(flow)

    xs = collect(range(xmin, xmax; length=length(axes(itp, 1))))
    ys = collect(range(ymin, ymax; length=length(axes(itp, 2))))

    hx_min = minimum(diff(xs))
    hy_min = minimum(diff(ys))
    domain_scale = max(
        abs(Float64(xmin)),
        abs(Float64(xmax)),
        abs(Float64(ymin)),
        abs(Float64(ymax)),
        1.0
    )

    boundary_tol = max(10 * tol, 100 * eps(Float64) * domain_scale)

    mask_boundary_filter_tol =
    min_mask_boundary_distance_cells * max(hx_min, hy_min)

    mask_segs_for_filter =
        ignore_masked_cells && min_mask_boundary_distance_cells > 0 ?
        maskboundarysegments(flow; zero_cell_tol=zero_cell_tol) :
        BoundarySegment{Float64,2}[]

    function distancetomaskboundary(p)
        dmin = Inf

        for seg in mask_segs_for_filter
            d = _pointsegmentdistance(p, seg.p0, seg.p1)
            dmin = min(dmin, d)
        end

        return dmin
    end

    function ismasked(v)
        return norm(v) < zero_cell_tol
    end

    candidates = Tuple{SVector{2,Float64}, SMatrix{2,2,Float64,4}}[]

    function quadraticroots(a, b, c; qtol=1e-12)
        roots = Float64[]

        if abs(a) < qtol
            abs(b) < qtol && return roots
            push!(roots, -c / b)
            return roots
        end

        Δ = b^2 - 4a*c
        Δ < -qtol && return roots
        Δ = max(Δ, 0.0)

        push!(roots, (-b + sqrt(Δ)) / (2a))
        push!(roots, (-b - sqrt(Δ)) / (2a))

        return roots
    end

    function tryaddroot!(ξ, η, x0, x1, y0, y1, a, b, v00, v10, v01, v11; local_tol=1e-9)
        if -local_tol <= ξ <= 1 + local_tol &&
           -local_tol <= η <= 1 + local_tol

            ξ = clamp(ξ, 0.0, 1.0)
            η = clamp(η, 0.0, 1.0)

            u = a[1] + a[2]*ξ + a[3]*η + a[4]*ξ*η
            v = b[1] + b[2]*ξ + b[3]*η + b[4]*ξ*η

            norm(SVector(u, v)) < max(1e-7, 10tol) || return

            x = x0 + ξ * (x1 - x0)
            y = y0 + η * (y1 - y0)
            p = SVector{2,Float64}(x, y)

            if ignore_masked_cells && min_mask_boundary_distance_cells > 0
                distancetomaskboundary(p) < mask_boundary_filter_tol && return
            end

            if ignore_boundary_points
                if p[1] <= xmin + boundary_tol ||
                   p[1] >= xmax - boundary_tol ||
                   p[2] <= ymin + boundary_tol ||
                   p[2] >= ymax - boundary_tol
                    return
                end
            end

            hx = x1 - x0
            hy = y1 - y0

            du_dξ = a[2] + a[4]*η
            du_dη = a[3] + a[4]*ξ
            dv_dξ = b[2] + b[4]*η
            dv_dη = b[3] + b[4]*ξ

            J = @SMatrix [
                du_dξ / hx  du_dη / hy
                dv_dξ / hx  dv_dη / hy
            ]

            duplicate = any(c -> norm(c[1] - p) < duplicate_tol, candidates)
            duplicate || push!(candidates, (p, J))
        end
    end

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

            if ignore_masked_cells
                nzero = count(v -> ismasked(v), (v00, v10, v01, v11))
                nzero >= 3 && continue
            end

            # u(ξ,η) = a1 + a2*ξ + a3*η + a4*ξ*η
            a = (
                v00[1],
                v10[1] - v00[1],
                v01[1] - v00[1],
                v11[1] - v10[1] - v01[1] + v00[1]
            )

            # v(ξ,η) = b1 + b2*ξ + b3*η + b4*ξ*η
            b = (
                v00[2],
                v10[2] - v00[2],
                v01[2] - v00[2],
                v11[2] - v10[2] - v01[2] + v00[2]
            )

            # Eliminate η:
            # η = -(a1 + a2*ξ) / (a3 + a4*ξ)
            q2 = b[2]*a[4] - b[4]*a[2]
            q1 = b[1]*a[4] + b[2]*a[3] - b[3]*a[2] - b[4]*a[1]
            q0 = b[1]*a[3] - b[3]*a[1]

            for ξ in quadraticroots(q2, q1, q0)
                denom = a[3] + a[4]*ξ
                abs(denom) < 1e-12 && continue

                η = -(a[1] + a[2]*ξ) / denom
                tryaddroot!(ξ, η, x0, x1, y0, y1, a, b, v00, v10, v01, v11)
            end

            # Fallback: eliminate ξ instead.
            r2 = b[3]*a[4] - b[4]*a[3]
            r1 = b[1]*a[4] + b[3]*a[2] - b[2]*a[3] - b[4]*a[1]
            r0 = b[1]*a[2] - b[2]*a[1]

            for η in quadraticroots(r2, r1, r0)
                denom = a[2] + a[4]*η
                abs(denom) < 1e-12 && continue

                ξ = -(a[1] + a[3]*η) / denom
                tryaddroot!(ξ, η, x0, x1, y0, y1, a, b, v00, v10, v01, v11)
            end
        end
    end

    cps = CriticalPoint{Float64,2}[]

    for (x, J) in candidates
        duplicate = any(cp -> norm(cp.x - x) < duplicate_tol, cps)
        duplicate && continue

        λ = eigvals(Matrix(J))
        kind = _classifyeigenvalues(λ)

        push!(cps, CriticalPoint(x, kind))
    end

    return cps
end

criticaltype(cp::CriticalPoint) = cp.kind

"""
    boundarybehavior(flow, t, x, normal)

Classify the boundary behavior at x using the outward normal.
Returns :inflow, :outflow, or :tangent.
"""
function boundarybehavior(flow::VCFlowData.InterpolatedFlow, x, normal; tol=1e-10)
    v = _flowvalue(flow, x)
    s = dot(v, normal)

    s < -tol && return :inflow
    s >  tol && return :outflow
    return :tangent
end

"""
    boundarysegments(flow; m=200)

Split the rectangular boundary into segments of uniform boundary behavior.
"""
function boundarysegments(flow::VCFlowData.InterpolatedFlow)
    itp = flow.itp
    xsamples = length(axes(itp, 1))
    ysamples = length(axes(itp, 2))
    xmin, ymin, xmax, ymax = _spatialbounds(flow)
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
        labels = [boundarybehavior(flow, p, normal) for p in pts]

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

function _ismaskedvalue(v; zero_cell_tol=1e-10)
    return norm(v) < zero_cell_tol
end


function maskboundarysegments(flow::VCFlowData.InterpolatedFlow; zero_cell_tol=1e-10)
    itp = flow.itp
    xmin, ymin, xmax, ymax = _spatialbounds(flow)

    nx = length(axes(itp, 1))
    ny = length(axes(itp, 2))

    xs = collect(range(Float64(xmin), Float64(xmax); length=nx))
    ys = collect(range(Float64(ymin), Float64(ymax); length=ny))

    segs = BoundarySegment{Float64,2}[]

    for i in 1:(nx - 1), j in 1:(ny - 1)
        x0 = xs[i]
        x1 = xs[i + 1]
        y0 = ys[j]
        y1 = ys[j + 1]

        xm = 0.5 * (x0 + x1)
        ym = 0.5 * (y0 + y1)

        # vertikale Interface-Prüfung: links <-> rechts
        p_left  = SVector{2,Float64}(x0, ym)
        p_right = SVector{2,Float64}(x1, ym)

        v_left  = _flowvalue(flow, p_left)
        v_right = _flowvalue(flow, p_right)

        m_left  = _ismaskedvalue(v_left; zero_cell_tol=zero_cell_tol)
        m_right = _ismaskedvalue(v_right; zero_cell_tol=zero_cell_tol)

        if m_left != m_right
            p0 = SVector{2,Float64}(xm, y0)
            p1 = SVector{2,Float64}(xm, y1)

            # Normal zeigt von Maske in Richtung Fluid
            normal = m_left ?
                SVector{2,Float64}(1.0, 0.0) :
                SVector{2,Float64}(-1.0, 0.0)

            push!(segs, BoundarySegment(p0, p1, normal))
        end

        # horizontale Interface-Prüfung: unten <-> oben
        p_bottom = SVector{2,Float64}(xm, y0)
        p_top    = SVector{2,Float64}(xm, y1)

        v_bottom = _flowvalue(flow, p_bottom)
        v_top    = _flowvalue(flow, p_top)

        m_bottom = _ismaskedvalue(v_bottom; zero_cell_tol=zero_cell_tol)
        m_top    = _ismaskedvalue(v_top; zero_cell_tol=zero_cell_tol)

        if m_bottom != m_top
            p0 = SVector{2,Float64}(x0, ym)
            p1 = SVector{2,Float64}(x1, ym)

            # Normal zeigt von Maske in Richtung Fluid
            normal = m_bottom ?
                SVector{2,Float64}(0.0, 1.0) :
                SVector{2,Float64}(0.0, -1.0)

            push!(segs, BoundarySegment(p0, p1, normal))
        end
    end

    return segs
end


function addswitchpointsfromsamples!(pts, xs, svals, normal, tangent; tol=1e-10)
    n = length(xs)
    n < 2 && return

    signs = map(s -> s > tol ? 1 : s < -tol ? -1 : 0, svals)

    i = 1
    while i < n
        # echter Vorzeichenwechsel zwischen zwei nicht-null Samples
        if signs[i] != 0 && signs[i + 1] != 0
            if signs[i] != signs[i + 1]
                s0 = svals[i]
                s1 = svals[i + 1]

                α = clamp(s0 / (s0 - s1), 0.0, 1.0)
                xsw = (1 - α) * xs[i] + α * xs[i + 1]

                push!(pts, BoundarySwitchPoint(xsw, normal, tangent))
            end

            i += 1
            continue
        end

        # zusammenhängende near-zero-Zone behandeln
        if signs[i] == 0
            run_start = i

            while i <= n && signs[i] == 0
                i += 1
            end

            run_end = i - 1

            left_sign = run_start > 1 ? signs[run_start - 1] : 0
            right_sign = i <= n ? signs[i] : 0

            # near-zero-Zone nur behalten, wenn links/rechts echte unterschiedliche Vorzeichen stehen
            if left_sign != 0 && right_sign != 0 && left_sign != right_sign
                mid = div(run_start + run_end, 2)
                push!(pts, BoundarySwitchPoint(xs[mid], normal, tangent))
            end

            continue
        end

        # signs[i] != 0 und signs[i+1] == 0
        i += 1
    end
end

"""
    boundaryswitchpoints(flow; m=400, tol=1e-10)

Find boundary switch points by detecting sign changes of v·n along each edge.
"""
function boundaryswitchpoints(
    flow::VCFlowData.InterpolatedFlow;
    tol=1e-10,
    patch::Bool=false,
    xsamples::Int=450,
    ysamples::Int=150,
    include_mask_boundary::Bool=true,
    zero_cell_tol=1e-10,
    duplicate_tol=1e-4,
    mask_sample_offset_cells::Real=0.25
)
    xmin, ymin, xmax, ymax = _spatialbounds(flow)

    T = Float64

    pts = BoundarySwitchPoint{Float64,2}[]

    # Äußere rechteckige Boundary

    if !patch
        edges = (
            (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmin, y), SVector{2,T}(-1, 0)), # left
            (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmax, y), SVector{2,T}( 1, 0)), # right
            (_linspace(xmin, xmax, xsamples), x -> SVector{2,T}(x, ymin), SVector{2,T}( 0,-1)), # bottom
            (_linspace(xmin, xmax, xsamples), x -> SVector{2,T}(x, ymax), SVector{2,T}( 0, 1)), # top
        )

        for (params, mkpt, normal) in edges
            tangent = SVector{2,T}(-normal[2], normal[1])

            xs = [mkpt(p) for p in params]
            svals = [dot(_flowvalue(flow, x), normal) for x in xs]

            addswitchpointsfromsamples!(pts, xs, svals, normal, tangent; tol=tol)
        end
    else
        side_edges = (
            (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmin, y), SVector{2,T}(-1, 0)), # left
            (_linspace(ymin, ymax, ysamples), y -> SVector{2,T}(xmax, y), SVector{2,T}( 1, 0)), # right
        )

        for (params, mkpt, normal) in side_edges
            tangent = SVector{2,T}(-normal[2], normal[1])

            xs = [mkpt(p) for p in params]
            svals = [dot(_flowvalue(flow, x), normal) for x in xs]

            addswitchpointsfromsamples!(pts, xs, svals, normal, tangent; tol=tol)
        end

        # top/bottom gepatcht über innere Sample-Linie und Tangentialkomponente
        xparams = _linspace(xmin, xmax, xsamples)
        yparams = _linspace(ymin, ymax, ysamples)

        if length(yparams) >= 3
            y_bottom_inner = yparams[2]
            y_top_inner = yparams[end - 1]

            patch_edges = (
                (
                    xparams,
                    x -> SVector{2,T}(x, y_bottom_inner),
                    x -> SVector{2,T}(x, ymin),
                    SVector{2,T}(0, -1)
                ),
                (
                    xparams,
                    x -> SVector{2,T}(x, y_top_inner),
                    x -> SVector{2,T}(x, ymax),
                    SVector{2,T}(0, 1)
                ),
            )

            for (params, offset_mkpt, boundary_mkpt, normal) in patch_edges
                tangent = SVector{2,T}(-normal[2], normal[1])

                xs_offset = [offset_mkpt(p) for p in params]
                xs_boundary = [boundary_mkpt(p) for p in params]

                svals = [dot(_flowvalue(flow, x), tangent) for x in xs_offset]

                addswitchpointsfromsamples!(pts, xs_boundary, svals, normal, tangent; tol=tol)
            end
        end
    end

    # Innere Maskengrenzen zusätzlich auswerten

    if include_mask_boundary
        mask_segs = maskboundarysegments(flow; zero_cell_tol=zero_cell_tol)

        itp = flow.itp
        nx = length(axes(itp, 1))
        ny = length(axes(itp, 2))

        xs_grid = collect(range(Float64(xmin), Float64(xmax); length=nx))
        ys_grid = collect(range(Float64(ymin), Float64(ymax); length=ny))

        hx = minimum(diff(xs_grid))
        hy = minimum(diff(ys_grid))
        h = min(hx, hy)

        sample_offset = mask_sample_offset_cells * h

        for seg in mask_segs
            normal = _safenormalize(seg.normal)
            tangent = _safenormalize(seg.p1 - seg.p0)

            params = range(0.0, 1.0; length=5)

            xs_boundary = [(1 - α) * seg.p0 + α * seg.p1 for α in params]

            # leicht in das Fluid hinein samplen, nicht direkt auf der Maskengrenze
            xs_probe = [x + sample_offset * normal for x in xs_boundary]

            svals = Float64[]

            for x in xs_probe
                if _inside(flow, x)
                    push!(svals, dot(_flowvalue(flow, x), normal))
                else
                    push!(svals, 0.0)
                end
            end

            addswitchpointsfromsamples!(pts, xs_boundary, svals, normal, tangent; tol=tol)
        end
    end


    # Deduplizieren
    out = BoundarySwitchPoint{Float64,2}[]

    for p in pts
        duplicate = any(q -> norm(q.x - p.x) < duplicate_tol, out)
        duplicate || push!(out, p)
    end

    return out
end

"""
    separatrixseeds(flow, cp::CriticalPoint; ϵ=1e-6)

For saddles:
- unstable directions -> :forward
- stable directions   -> :backward
"""
function separatrixseeds(flow::VCFlowData.InterpolatedFlow, cp::CriticalPoint; ϵ=1e-6)
    cp.kind isa Saddle || return Tuple{SVector{2,Float64},Symbol}[]

    x0 = SVector{2,Float64}(cp.x)
    J = ForwardDiff.jacobian(x0 -> _flowvalue(flow, x0), x0)
    F = eigen(Matrix(J))

    seeds = Tuple{SVector{2,Float64},Symbol}[]

    for i in eachindex(F.values)
        λ = F.values[i]
        v = SVector{2,Float64}(F.vectors[:, i])
        v = _safenormalize(v)
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
    separatrixseeds(flow, seg::BoundarySegment; k=20, ϵ=1e-6)

Sample inflow boundary segments and move seeds slightly into the domain.
"""
function separatrixseeds(flow::VCFlowData.InterpolatedFlow, seg::BoundarySegment; k=20, ϵ=1e-6)
    mid = 0.5 * (seg.p0 + seg.p1)
    beh = boundarybehavior(flow, mid, seg.normal)

    beh == :inflow || return Tuple{SVector{2,Float64},Symbol}[]

    n̂ = _safenormalize(SVector{2,Float64}(seg.normal))
    seeds = Tuple{SVector{2,Float64},Symbol}[]

    for s in range(0.0, 1.0; length=k)
        x = (1-s) * SVector{2,Float64}(seg.p0) + s * SVector{2,Float64}(seg.p1)
        x_in = x - ϵ * n̂
        push!(seeds, (x_in, :forward))
    end

    return seeds
end

"""
    separatrixseeds(flow, bsp::BoundarySwitchPoint; ϵ=1e-3)

Create one seed slightly inside the domain near a boundary switch point.
From this seed, the separatrix should be traced both forward and backward.
"""
function separatrixseeds(
    flow::VCFlowData.InterpolatedFlow,
    bsp::BoundarySwitchPoint;
    ϵ=1e-3
)
    x0 = SVector{2,Float64}(bsp.x)
    n̂ = _safenormalize(SVector{2,Float64}(bsp.normal))

    # move only slightly into the domain
    xin = x0 - ϵ * n̂

    return [
        (xin, :forward),
        (xin, :backward),
    ]
end

integrationdirection(::BoundarySegment) = :forward

"""
    traceseparatrix(flow, x0; dir=:forward, h=0.005, maxsteps=4000,
                     stop_eps=5e-3, minsteps_before_stop=10)

Trace a separatrix as sampled pathline using RK43.
"""
function traceseparatrix(flow::VCFlowData.InterpolatedFlow, x0::SVector{2,Float64};
    dir::Symbol=:forward,
    h::Float64=0.005,
    maxsteps::Int=4000,
    stop_eps::Float64=5e-3,
    minsteps_before_stop::Int=10
)
    cps = criticalpoints(flow)
    bsps = boundaryswitchpoints(flow)

    s = dir === :forward ? 1.0 : -1.0

    function dy(t, y)
        _inside(flow, y) || return RK43.OutOfDomain
        return s * _flowvalue(flow, y)
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
    tracesaddleseparatrix(flow, cp, x0; dir=:forward, h=0.005,
                            maxsteps=4000, stop_eps=1e-2, minsteps_before_stop=10)

Trace a separatrix starting from a saddle seed.

Delegates the actual integration to `traceseparatrix`.
"""
function tracesaddleseparatrix(
    flow::VCFlowData.InterpolatedFlow,
    cp::CriticalPoint,
    x0::SVector{2,Float64};
    dir::Symbol=:forward,
    h::Float64=0.005,
    maxsteps::Int=4000,
    stop_eps::Float64=1e-2,
    minsteps_before_stop::Int=10
)
    return traceseparatrix(flow, x0;
        dir=dir,
        h=h,
        maxsteps=maxsteps,
        stop_eps=stop_eps,
        minsteps_before_stop=minsteps_before_stop
    )
end

function integrationdirection(cp::CriticalPoint)
    cp.kind isa Source && return :backward
    cp.kind isa Sink   && return :forward
    cp.kind isa Saddle && return :both
    return :forward
end

export criticalpoints, criticaltype, boundarysegments, boundaryswitchpoints, separatrixseeds, divergence