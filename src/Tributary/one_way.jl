# =============================================================================
# One-Way Directed Tributary Areas
# =============================================================================

"""
    get_tributary_polygons_one_way(vertices; weights=nothing, axis, buffers=nothing, holes=nothing)

Compute tributary polygons using one-way directed partitioning along the specified axis.

Load travels parallel to `axis`. The floor is swept in strips across `axis`;
inside a strip each stretch of floor runs between two boundary crossings and
goes to the edges at its ends, split by `weights`.

## Arguments
- `vertices::Vector{<:Point}`: Polygon vertices as Meshes.Point objects
- `weights::Union{Nothing, AbstractVector{<:Real}}`: Optional edge weights (one per edge)
- `axis::AbstractVector{<:Real}`: Direction vector [vx, vy] for load distribution
- `buffers`: Accepted for compatibility and not used.
- `holes`: Optional openings, each a vertex loop (`Point`s or `(x, y)` in meters).
  An opening edge is not a support: a stretch that ends on an opening goes
  wholly to the edge at its other end. A stretch with an opening at both ends
  has no support along `axis`, and this throws.

## Returns
`Vector{TributaryPolygon}` in parametric form. Edges parallel to axis get zero area.
Use `vertices(trib, beam_start, beam_end)` to get absolute coordinates.

Each profile is measured along `axis`: a stretch of depth `D` (along `axis`)
reaching an edge whose direction makes angle θ with `axis` loads it with
`D·|sin θ|` per unit edge length. That keeps the whole stretch on the edge it
reaches, including skewed edges. Areas are net of openings and sum to the net
floor area; this throws otherwise.
"""
function get_tributary_polygons_one_way(
    vertices::Vector{<:Point};
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    axis::AbstractVector{<:Real},
    buffers::Union{Nothing, TributaryBuffers} = nothing,
    holes = nothing
)
    m = length(vertices)
    m >= 3 || return TributaryPolygon[]
    
    # Convert to 2D coords and ensure CCW
    pts_orig = [_to_2d(v) for v in vertices]
    pts_orig = _ensure_ccw(pts_orig)
    _, hole_loops = _span_region(pts_orig, holes)
    
    # Normalize axis: u = load direction, n = perpendicular (sweep direction)
    vx, vy = Float64(axis[1]), Float64(axis[2])
    vlen = hypot(vx, vy)
    vlen < 1e-12 && error("axis must be non-zero")
    u = (vx / vlen, vy / vlen)
    n = (-u[2], u[1])
    
    # Handle weights
    w = isnothing(weights) ? ones(m) : Float64.(weights)
    length(w) == m || error("weights must have same length as vertices ($m)")
    
    # Transform to (s,t) coordinates
    pts_st = [_to_st(p, u, n) for p in pts_orig]
    holes_st = [NTuple{2,Float64}[_to_st(p, u, n) for p in h] for h in hole_loops]

    # Boundary segments in (s, t). Outer edge i has id i; opening edges have id 0.
    segs = Tuple{NTuple{2,Float64}, NTuple{2,Float64}, Int}[]
    for i in 1:m
        push!(segs, (pts_st[i], pts_st[mod1(i + 1, m)], i))
    end
    for h in holes_st, j in eachindex(h)
        push!(segs, (h[j], h[mod1(j + 1, length(h))], 0))
    end
    
    # Strip boundaries: every vertex, and every place an opening crosses the
    # outline or another opening, so each strip has one arrangement of stretches.
    swap(loop) = NTuple{2,Float64}[(p[2], p[1]) for p in loop]
    t_lo = minimum(p[2] for p in pts_st)
    t_hi = maximum(p[2] for p in pts_st)
    t_all = Float64[p[2] for p in pts_st]
    for (k, h) in enumerate(holes_st)
        append!(t_all, (p[2] for p in h))
        append!(t_all, _crossing_xs(swap(pts_st), swap(h)))
        for h2 in holes_st[k+1:end]
            append!(t_all, _crossing_xs(swap(h), swap(h2)))
        end
    end
    t_vals = sort!(unique!(clamp.(t_all, t_lo, t_hi)))
    empty_result() = [_make_tributary_with_depths(i, NTuple{2,Float64}[], pts_orig, 0.0, 0.0) for i in 1:m]
    length(t_vals) < 2 && return empty_result()
    
    edge_areas = zeros(m)
    edge_pieces = [NTuple{4,Float64}[] for _ in 1:m]     # (α0, w0, α1, w1) along the edge
    edge_interior_pts = [NTuple{2,Float64}[] for _ in 1:m]
    edge_t_range = [(Inf, -Inf) for _ in 1:m]
    edge_is_left = [true for _ in 1:m]

    function credit!(e, d0, d1, t0, t1, interior0, interior1, is_left)
        (d0 > 1e-12 || d1 > 1e-12) || return
        edge_areas[e] += (d0 + d1) / 2 * (t1 - t0)
        (sa, ta) = pts_st[e]
        (sb, tb) = pts_st[mod1(e + 1, m)]
        Δt = tb - ta
        k = abs(Δt) / hypot(sb - sa, Δt)
        push!(edge_pieces[e], ((t0 - ta) / Δt, d0 * k, (t1 - ta) / Δt, d1 * k))
        push!(edge_interior_pts[e], interior0)
        push!(edge_interior_pts[e], interior1)
        edge_is_left[e] = is_left
        tmin, tmax = edge_t_range[e]
        edge_t_range[e] = (min(tmin, t0), max(tmax, t1))
    end
    
    for k in 1:(length(t_vals) - 1)
        t0, t1 = t_vals[k], t_vals[k + 1]
        t1 - t0 < 1e-12 && continue
        for (iL, iR) in _strip_stretches(segs, (t0 + t1) / 2)
            sL0, sL1 = _seg_s_at_t(segs[iL], t0), _seg_s_at_t(segs[iL], t1)
            sR0, sR1 = _seg_s_at_t(segs[iR], t0), _seg_s_at_t(segs[iR], t1)
            (sR0 - sL0 > 1e-9 || sR1 - sL1 > 1e-9) || continue
            eL, eR = segs[iL][3], segs[iR][3]
            if eL > 0 && eR > 0
                split0 = _weighted_split(sL0, sR0, w[eL], w[eR])
                split1 = _weighted_split(sL1, sR1, w[eL], w[eR])
            elseif eL > 0
                split0, split1 = sR0, sR1
            elseif eR > 0
                split0, split1 = sL0, sL1
            else
                throw(ArgumentError(
                    "One-way strip between t = $(round(t0, digits=4)) and $(round(t1, digits=4)) m " *
                    "runs from opening to opening along the span axis, so no edge supports it."))
            end
            eL > 0 && credit!(eL, max(0.0, split0 - sL0), max(0.0, split1 - sL1), t0, t1,
                              (split0, t0), (split1, t1), true)
            eR > 0 && credit!(eR, max(0.0, sR0 - split0), max(0.0, sR1 - split1), t0, t1,
                              (split0, t0), (split1, t1), false)
        end
    end
    
    net_area = _net_area(pts_orig, hole_loops)
    results = TributaryPolygon[]
    
    for i in 1:m
        if edge_areas[i] < 1e-12
            push!(results, _make_tributary_with_depths(i, NTuple{2,Float64}[], pts_orig, 0.0, 0.0))
            continue
        end
        
        # Get the actual polygon edge vertices for this edge
        v1_st = pts_st[i]
        v2_st = pts_st[mod1(i + 1, m)]
        
        # Clip the actual polygon edge to the contributing t-range
        t_min, t_max = edge_t_range[i]
        actual_exterior = _get_edge_in_t_range(v1_st, v2_st, t_min, t_max)
        
        # Reconstruct polygon from actual exterior + interior (split) points
        poly_xy = _build_tributary_polygon_v2(
            actual_exterior, edge_interior_pts[i], edge_is_left[i], u, n
        )
        
        area = edge_areas[i]
        frac = net_area > 0 ? area / net_area : 0.0
        push!(results, _make_tributary_with_depths(i, poly_xy, pts_orig, area, frac;
                                                   profile = _assemble_profile(edge_pieces[i])))
    end
    
    assert_tributary_balance(results, pts_orig, net_area; context = " (one-way)")
    return results
end

"""`s` where segment `seg = (p1, p2, id)` meets the line at `t` (clamped to the segment)."""
function _seg_s_at_t(seg, t::Float64)
    (s1, t1), (s2, t2), _ = seg
    abs(t2 - t1) < 1e-12 && return (s1 + s2) / 2
    α = clamp((t - t1) / (t2 - t1), 0.0, 1.0)
    return s1 + α * (s2 - s1)
end

"""
Stretches of floor on the line at `t`, as `(left, right)` indices into `segs`.
The outline's spans minus the openings' spans; openings do not overlap each other.
"""
function _strip_stretches(segs, t::Float64)
    outer = Tuple{Float64, Int}[]
    voids = Tuple{Float64, Int}[]
    for (k, ((s1, t1), (s2, t2), id)) in enumerate(segs)
        abs(t2 - t1) < 1e-12 && continue
        if (t1 <= t < t2) || (t2 <= t < t1)
            push!(id > 0 ? outer : voids, (s1 + (t - t1) / (t2 - t1) * (s2 - s1), k))
        end
    end
    sort!(outer; by = first)
    sort!(voids; by = first)
    stretches = Tuple{Int, Int}[]
    for a in 1:2:(length(outer) - 1)
        lo, hi = outer[a], outer[a + 1]
        for b in 1:2:(length(voids) - 1)
            c, d = voids[b], voids[b + 1]
            (d[1] <= lo[1] || c[1] >= hi[1]) && continue
            c[1] > lo[1] && push!(stretches, (lo[2], c[2]))
            lo = d
        end
        hi[1] > lo[1] && push!(stretches, (lo[2], hi[2]))
    end
    return stretches
end

"""
Width profile on `[0, 1]` from `(α0, w0, α1, w1)` pieces along one edge.
Pieces from different strips do not overlap. Gaps carry zero.
"""
function _assemble_profile(pieces::Vector{NTuple{4,Float64}})
    spans = [a0 <= a1 ? (a0, w0, a1, w1) : (a1, w1, a0, w0) for (a0, w0, a1, w1) in pieces]
    sort!(spans; by = first)
    pos = Float64[]
    w = Float64[]
    cursor = 0.0
    for (a, wa, b, wb) in spans
        a = max(clamp(a, 0.0, 1.0), cursor)
        b = clamp(b, 0.0, 1.0)
        b < a && continue
        if a > cursor + 1e-12
            push!(pos, cursor, a)
            push!(w, 0.0, 0.0)
        end
        if !isempty(pos) && pos[end] == a && w[end] == wa
            push!(pos, b)
            push!(w, wb)
        else
            push!(pos, a, b)
            push!(w, wa, wb)
        end
        cursor = b
    end
    if cursor < 1.0 - 1e-12
        push!(pos, cursor, 1.0)
        push!(w, 0.0, 0.0)
    end
    isempty(pos) && return ([0.0, 1.0], [0.0, 0.0])
    pos[1] = 0.0
    pos[end] = 1.0
    return pos, w
end

"""
Get points along a polygon edge within a t-range, including the actual vertices.
Returns points ordered by increasing t.
"""
function _get_edge_in_t_range(v1::NTuple{2,Float64}, v2::NTuple{2,Float64}, 
                               t_min::Float64, t_max::Float64)
    s1, t1 = v1
    s2, t2 = v2
    
    # Handle horizontal edge
    if abs(t2 - t1) < 1e-12
        if t_min - 1e-9 <= t1 <= t_max + 1e-9
            return t1 <= t2 ? [v1, v2] : [v2, v1]
        else
            return NTuple{2,Float64}[]
        end
    end
    
    # Collect points on edge within [t_min, t_max]
    points = NTuple{2,Float64}[]
    
    # Add actual vertices if within range
    if t_min - 1e-9 <= t1 <= t_max + 1e-9
        push!(points, v1)
    end
    if t_min - 1e-9 <= t2 <= t_max + 1e-9
        push!(points, v2)
    end
    
    # Add clipped endpoints at t_min and t_max if they're strictly inside edge's t-range
    edge_t_lo, edge_t_hi = minmax(t1, t2)
    
    if t_min > edge_t_lo + 1e-9 && t_min < edge_t_hi - 1e-9
        # t_min is strictly inside edge - add point at t_min
        α = (t_min - t1) / (t2 - t1)
        s_at_tmin = s1 + α * (s2 - s1)
        push!(points, (s_at_tmin, t_min))
    end
    
    if t_max > edge_t_lo + 1e-9 && t_max < edge_t_hi - 1e-9
        # t_max is strictly inside edge - add point at t_max
        α = (t_max - t1) / (t2 - t1)
        s_at_tmax = s1 + α * (s2 - s1)
        push!(points, (s_at_tmax, t_max))
    end
    
    # Sort by t
    sort!(points, by=p->p[2])
    
    return points
end

# =============================================================================
# Polygon Reconstruction with Exact Exterior Boundary
# =============================================================================

"""
Build tributary polygon from exterior (polygon edge) and interior (split) points.

For edge on LEFT of interval:
- Exterior is on LEFT, Interior is on RIGHT
- CCW order: exterior going UP, then interior going DOWN

For edge on RIGHT of interval:  
- Exterior is on RIGHT, Interior is on LEFT
- CCW order: interior going UP, then exterior going DOWN
"""
function _build_tributary_polygon_v2(exterior_pts::Vector{NTuple{2,Float64}},
                                     interior_pts::Vector{NTuple{2,Float64}},
                                     is_left_edge::Bool,
                                     u::NTuple{2,Float64}, n::NTuple{2,Float64})
    
    (isempty(exterior_pts) || isempty(interior_pts)) && return NTuple{2,Float64}[]
    
    # Sort and deduplicate both chains by t
    exterior_chain = _sort_and_dedup(exterior_pts)
    interior_chain = _sort_and_dedup(interior_pts)
    
    (length(exterior_chain) < 2 || length(interior_chain) < 2) && return NTuple{2,Float64}[]
    
    # Build polygon in CCW order
    poly_st = NTuple{2,Float64}[]
    
    if is_left_edge
        # Exterior on LEFT, Interior on RIGHT
        # CCW: start at bottom-left (exterior), go UP along exterior,
        #      then DOWN along interior (reversed), back to start
        
        # Exterior going up (bottom to top)
        for pt in exterior_chain
            if isempty(poly_st) || !_pts_equal(pt, poly_st[end])
                push!(poly_st, pt)
            end
        end
        
        # Interior going down (top to bottom)
        for i in length(interior_chain):-1:1
            pt = interior_chain[i]
            if isempty(poly_st) || !_pts_equal(pt, poly_st[end])
                push!(poly_st, pt)
            end
        end
    else
        # Exterior on RIGHT, Interior on LEFT
        # CCW: start at bottom-left (interior), go UP along interior,
        #      then DOWN along exterior (reversed), back to start
        
        # Interior going up (bottom to top)
        for pt in interior_chain
            if isempty(poly_st) || !_pts_equal(pt, poly_st[end])
                push!(poly_st, pt)
            end
        end
        
        # Exterior going down (top to bottom)
        for i in length(exterior_chain):-1:1
            pt = exterior_chain[i]
            if isempty(poly_st) || !_pts_equal(pt, poly_st[end])
                push!(poly_st, pt)
            end
        end
    end
    
    length(poly_st) < 3 && return NTuple{2,Float64}[]
    
    # Close polygon
    if !_pts_equal(poly_st[1], poly_st[end])
        push!(poly_st, poly_st[1])
    end
    
    # Transform to (x,y)
    poly_xy = [_to_xy(s, t, u, n) for (s, t) in poly_st]
    
    # Simplify and ensure CCW
    poly_xy = _simplify_collinear(poly_xy)
    poly_xy = _ensure_ccw(poly_xy)
    
    return poly_xy
end

"""Sort points by t and remove duplicates."""
function _sort_and_dedup(pts::Vector{NTuple{2,Float64}})
    sorted = sort(pts, by=p->p[2])
    result = NTuple{2,Float64}[]
    for pt in sorted
        if isempty(result) || !_pts_equal(pt, result[end])
            push!(result, pt)
        end
    end
    return result
end

# =============================================================================
# Coordinate Transforms
# =============================================================================

_to_st(p::NTuple{2,Float64}, u::NTuple{2,Float64}, n::NTuple{2,Float64}) = 
    (p[1]*u[1] + p[2]*u[2], p[1]*n[1] + p[2]*n[2])

_to_xy(s::Float64, t::Float64, u::NTuple{2,Float64}, n::NTuple{2,Float64}) = 
    (s*u[1] + t*n[1], s*u[2] + t*n[2])

# =============================================================================
# Weighted Split
# =============================================================================

"""Compute weighted split point between sL and sR based on edge weights."""
function _weighted_split(sL::Float64, sR::Float64, wL::Float64, wR::Float64)
    denom = wL + wR
    denom < 1e-12 && return (sL + sR) / 2
    α = wR / denom
    return (1 - α) * sL + α * sR
end

# =============================================================================
# Polygon Utilities
# =============================================================================

"""Remove collinear points from polygon."""
function _simplify_collinear(pts::Vector{NTuple{2,Float64}})
    length(pts) < 4 && return pts
    
    closed = _pts_equal(pts[1], pts[end])
    work = closed ? pts[1:end-1] : pts
    
    nv = length(work)
    nv < 3 && return pts
    
    result = NTuple{2,Float64}[]
    
    for i in 1:nv
        prev = work[mod1(i - 1, nv)]
        curr = work[i]
        next = work[mod1(i + 1, nv)]
        
        cross = (curr[1] - prev[1]) * (next[2] - curr[2]) - 
                (curr[2] - prev[2]) * (next[1] - curr[1])
        
        if abs(cross) > 1e-9
            push!(result, curr)
        end
    end
    
    length(result) < 3 && return pts
    
    push!(result, result[1])
    return result
end

_pts_equal(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = hypot(a[1]-b[1], a[2]-b[2]) < 1e-9
