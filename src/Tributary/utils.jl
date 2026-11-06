# Tributary Area Types and Utilities

import Meshes: Point, coords
using Unitful: ustrip, @u_str

# =============================================================================
# Reusable Buffer for Tributary Computation (reduces GC pressure)
# =============================================================================

"""
    TributaryBuffers

Pre-allocated buffers for tributary area computations.
Reuse across multiple cells to reduce allocations.

# Example
```julia
buffers = TributaryBuffers(max_vertices=20)
for cell in cells
    results = get_tributary_polygons(verts; buffers=buffers)
end
```
"""
mutable struct TributaryBuffers
    # 2D coordinate buffers
    pts_2d::Vector{NTuple{2,Float64}}
    pts_st::Vector{NTuple{2,Float64}}
    # Weight buffers
    weights::Vector{Float64}
    # Edge tracking buffers  
    edge_areas::Vector{Float64}
    edge_interior_pts::Vector{Vector{NTuple{2,Float64}}}
    edge_t_range::Vector{Tuple{Float64,Float64}}
    edge_is_left::Vector{Bool}
    # Scanline working buffers
    crosses::Vector{Tuple{Float64,Int}}
    intervals::Vector{Tuple{Tuple{Float64,Int}, Tuple{Float64,Int}}}
    # Polygon reconstruction
    poly_st::Vector{NTuple{2,Float64}}
    poly_xy::Vector{NTuple{2,Float64}}
end

"""Create buffers sized for polygons with up to `max_vertices` vertices."""
function TributaryBuffers(; max_vertices::Int=32)
    TributaryBuffers(
        Vector{NTuple{2,Float64}}(undef, max_vertices),
        Vector{NTuple{2,Float64}}(undef, max_vertices),
        Vector{Float64}(undef, max_vertices),
        Vector{Float64}(undef, max_vertices),
        [NTuple{2,Float64}[] for _ in 1:max_vertices],
        Vector{Tuple{Float64,Float64}}(undef, max_vertices),
        Vector{Bool}(undef, max_vertices),
        Tuple{Float64,Int}[],
        Tuple{Tuple{Float64,Int}, Tuple{Float64,Int}}[],
        NTuple{2,Float64}[],
        NTuple{2,Float64}[]
    )
end

"""Ensure buffers can handle `n` vertices, growing if needed."""
function ensure_capacity!(buf::TributaryBuffers, n::Int)
    if length(buf.pts_2d) < n
        resize!(buf.pts_2d, n)
        resize!(buf.pts_st, n)
        resize!(buf.weights, n)
        resize!(buf.edge_areas, n)
        while length(buf.edge_interior_pts) < n
            push!(buf.edge_interior_pts, NTuple{2,Float64}[])
        end
        resize!(buf.edge_t_range, n)
        resize!(buf.edge_is_left, n)
    end
    # Reset accumulators
    fill!(buf.edge_areas, 0.0)
    for v in buf.edge_interior_pts
        empty!(v)
    end
    fill!(buf.edge_t_range, (Inf, -Inf))
    fill!(buf.edge_is_left, true)
    empty!(buf.crosses)
    empty!(buf.intervals)
    empty!(buf.poly_st)
    empty!(buf.poly_xy)
end

"""
Parametric tributary polygon relative to a beam edge.

All length values are in **meters** (SI base unit).

Fields:
- `local_edge_idx`: Edge index within cell (1..n_edges)
- `s`: Normalized positions along beam [0,1] (unitless)
- `d`: Perpendicular distances from beam (meters) - these are TRIBUTARY depths
- `area`: Tributary area (m²)
- `fraction`: Fraction of total cell area (unitless)

Cell depth fields (perpendicular distance to CELL boundary, for EFM l2 calculation):
- `cell_depths_s`: Parametric positions along edge where depths were sampled
- `cell_depths`: Perpendicular depths at each sample point (m)
- `l2_stiff`: Stiffness-consistent effective width via cubic mean (m) - PRIMARY for EFM
- `cell_depth_max`: Maximum perpendicular depth (quick accessor for max span)
- `cell_depth_s_max`: Parametric position where max depth occurs

Load profile (what the edge carries):
- `profile_s`: Positions along the edge, nondecreasing from 0 to 1. A repeated
  position is a step.
- `profile_w`: Tributary width at each position (m), piecewise linear between
  positions. `∫ w ds × edge length = area` exactly.

The partition writes the profile. The `(s, d)` outline is for drawing, flange
widths, and strips; it is not a load profile. A one-way strip reaches a skewed
edge at an angle, so its outline is not a width over the edge.
"""
struct TributaryPolygon
    local_edge_idx::Int
    s::Vector{Float64}
    d::Vector{Float64}
    area::Float64
    fraction::Float64
    # Cell boundary depths (for EFM l2 calculation)
    cell_depths_s::Vector{Float64}  # parametric positions of samples
    cell_depths::Vector{Float64}    # depths at each sample (m)
    l2_stiff::Float64               # cubic mean effective width (m)
    cell_depth_max::Float64         # max depth for quick access
    cell_depth_s_max::Float64       # s position of max
    profile_s::Vector{Float64}
    profile_w::Vector{Float64}
end

"""
Outline-only constructor. The profile is the outline sliced across the edge,
which is the load profile when every slice starts on the edge: two-way faces,
and hand-built test polygons.
"""
function TributaryPolygon(local_edge_idx::Int, s::Vector{Float64}, d::Vector{Float64},
                          area::Float64, fraction::Float64,
                          cell_depths_s::Vector{Float64}, cell_depths::Vector{Float64},
                          l2_stiff::Float64, cell_depth_max::Float64, cell_depth_s_max::Float64)
    ps, pw = _outline_profile(s, d)
    TributaryPolygon(local_edge_idx, s, d, area, fraction,
                     cell_depths_s, cell_depths, l2_stiff, cell_depth_max, cell_depth_s_max,
                     ps, pw)
end

# Constructor with default cell depths (for backwards compatibility)
function TributaryPolygon(local_edge_idx::Int, s::Vector{Float64}, d::Vector{Float64}, 
                          area::Float64, fraction::Float64)
    TributaryPolygon(local_edge_idx, s, d, area, fraction, 
                     Float64[], Float64[], 0.0, 0.0, 0.0)
end

"""The same tributary with `fraction` replaced."""
_with_fraction(t::TributaryPolygon, f::Float64) =
    TributaryPolygon(t.local_edge_idx, t.s, t.d, t.area, f,
                     t.cell_depths_s, t.cell_depths, t.l2_stiff, t.cell_depth_max,
                     t.cell_depth_s_max, t.profile_s, t.profile_w)

"""Area under a tributary's load profile for an edge of length `edge_len` (m²)."""
function profile_area(trib::TributaryPolygon, edge_len::Real)
    return _trapezoid(trib.profile_s, trib.profile_w) * edge_len
end

"""`∫ w ds` over a piecewise-linear profile."""
function _trapezoid(pos::Vector{Float64}, w::Vector{Float64})
    acc = 0.0
    @inbounds for i in 1:(length(pos) - 1)
        acc += (w[i] + w[i + 1]) * 0.5 * (pos[i + 1] - pos[i])
    end
    return acc
end

"""Profile of an `(s, d)` outline: its slice length across the edge, on `s ∈ [0, 1]`."""
function _outline_profile(s::Vector{Float64}, d::Vector{Float64})
    isempty(s) && return ([0.0, 1.0], [0.0, 0.0])
    poly = NTuple{2,Float64}[(s[i], d[i]) for i in eachindex(s)]
    pos, w, _ = _edge_slice_profile(poly, Vector{NTuple{2,Float64}}[])
    return pos, w
end

# =============================================================================
# Slicing: exact width of a region along lines x = const
# =============================================================================

"""Intervals of the line `x = x0` inside the polygon `poly`, sorted (lo, hi)."""
function _vertical_intervals(poly::Vector{NTuple{2,Float64}}, x0::Float64)
    ys = Float64[]
    n = length(poly)
    @inbounds for i in 1:n
        x1, y1 = poly[i]
        x2, y2 = poly[mod1(i + 1, n)]
        if (x1 <= x0 < x2) || (x2 <= x0 < x1)
            push!(ys, y1 + (x0 - x1) * (y2 - y1) / (x2 - x1))
        end
    end
    sort!(ys)
    return [(ys[k], ys[k + 1]) for k in 1:2:(length(ys) - 1)]
end

"""Length of `x = x0` inside `poly` and outside every hole. Holes do not overlap each other."""
function _slice_length(poly, holes, x0::Float64)
    outer = _vertical_intervals(poly, x0)
    total = 0.0
    for (a, b) in outer
        total += b - a
    end
    for h in holes, (c, d) in _vertical_intervals(h, x0), (a, b) in outer
        total -= max(0.0, min(b, d) - max(a, c))
    end
    return max(total, 0.0)
end

"""x of every proper crossing between an edge of `p` and an edge of `q`."""
function _crossing_xs(p::Vector{NTuple{2,Float64}}, q::Vector{NTuple{2,Float64}})
    xs = Float64[]
    np, nq = length(p), length(q)
    for i in 1:np
        a = p[i]
        b = p[mod1(i + 1, np)]
        r = (b[1] - a[1], b[2] - a[2])
        for j in 1:nq
            c = q[j]
            d = q[mod1(j + 1, nq)]
            s = (d[1] - c[1], d[2] - c[2])
            den = r[1] * s[2] - r[2] * s[1]
            abs(den) < 1e-15 && continue
            w = (c[1] - a[1], c[2] - a[2])
            t = (w[1] * s[2] - w[2] * s[1]) / den
            u = (w[1] * r[2] - w[2] * r[1]) / den
            (0.0 <= t <= 1.0 && 0.0 <= u <= 1.0) && push!(xs, a[1] + t * r[1])
        end
    end
    return xs
end

"""
    _slice_profile(poly, holes, lo, hi) -> (positions, widths, inside, total, outside)

Width of `poly` minus `holes` along each line `x = const`, as an exact
piecewise-linear profile on `[lo, hi]`.

Between consecutive vertex and crossing abscissae the same edges cross every
slice, so the width is linear there. Each piece is sampled at its thirds and
extended to its ends, which gives the one-sided values at a step. `inside` is
`∫ width dx` over `[lo, hi]`; `total` is over the whole region. `outside`
holds the `(a, wa, b, wb)` pieces beyond `[lo, hi]`.
"""
function _slice_profile(poly::Vector{NTuple{2,Float64}},
                        holes::Vector{Vector{NTuple{2,Float64}}},
                        lo::Float64, hi::Float64)
    xs = Float64[lo, hi]
    for p in poly
        push!(xs, p[1])
    end
    for h in holes
        for p in h
            push!(xs, p[1])
        end
        append!(xs, _crossing_xs(poly, h))
    end
    sort!(xs)
    pos = Float64[]
    w = Float64[]
    outside = NTuple{4,Float64}[]
    inside = 0.0
    total = 0.0
    for k in 1:(length(xs) - 1)
        a, b = xs[k], xs[k + 1]
        Δ = b - a
        Δ < 1e-12 && continue
        f1 = _slice_length(poly, holes, a + Δ / 3)
        f2 = _slice_length(poly, holes, a + 2Δ / 3)
        fa = max(0.0, 2 * f1 - f2)
        fb = max(0.0, 2 * f2 - f1)
        piece = 0.5 * (fa + fb) * Δ
        total += piece
        if !(a >= lo - 1e-12 && b <= hi + 1e-12)
            piece > 0 && push!(outside, (a, fa, b, fb))
            continue
        end
        inside += piece
        if !isempty(pos) && pos[end] == a && w[end] == fa
            push!(pos, b)
            push!(w, fb)
        else
            push!(pos, a, b)
            push!(w, fa, fb)
        end
    end
    if isempty(pos)
        return ([lo, hi], [0.0, 0.0], 0.0, total, outside)
    end
    pos[1] > lo && (pushfirst!(pos, lo, pos[1]); pushfirst!(w, 0.0, 0.0))
    pos[end] < hi && (push!(pos, pos[end], hi); push!(w, 0.0, 0.0))
    pos[1] = lo
    pos[end] = hi
    return pos, w, inside, total, outside
end

"""
    _edge_slice_profile(poly, holes) -> (positions, widths, area)

[`_slice_profile`](@ref) of a region in its edge's frame, on `[0, 1]`. A part
past an end is mirrored back onto the edge about that end: a straight-skeleton
region runs past its edge at a re-entrant corner, and that load belongs at the
corner, so it stays on this edge beside the support. `area` is the whole
region's `∫ width ds`.
"""
function _edge_slice_profile(poly::Vector{NTuple{2,Float64}},
                             holes::Vector{Vector{NTuple{2,Float64}}})
    pos, w, _, total, outside = _slice_profile(poly, holes, 0.0, 1.0)
    isempty(outside) && return pos, w, total
    profiles = [(pos, w)]
    for (a, fa, b, fb) in outside
        a2, fa2, b2, fb2 = b <= 0.0 ? (-b, fb, -a, fa) : (2.0 - b, fb, 2.0 - a, fa)
        a2, b2 = clamp(a2, 0.0, 1.0), clamp(b2, 0.0, 1.0)
        push!(profiles, ([0.0, a2, a2, b2, b2, 1.0], [0.0, 0.0, fa2, fb2, 0.0, 0.0]))
    end
    p, ww = _sum_profiles(profiles)
    return p, ww, total
end

"""
    _sum_profiles(profiles) -> (positions, widths)

Sum piecewise-linear profiles on the union of their positions. A repeated
position is a step; the sum keeps both sides of every step.
"""
function _sum_profiles(profiles::Vector{Tuple{Vector{Float64}, Vector{Float64}}})
    length(profiles) == 1 && return profiles[1]
    breaks = sort!(unique!(reduce(vcat, (p for (p, _) in profiles))))
    pos = Float64[]
    w = Float64[]
    for s in breaks
        l = sum(_profile_limit(p, ww, s, true) for (p, ww) in profiles)
        r = sum(_profile_limit(p, ww, s, false) for (p, ww) in profiles)
        push!(pos, s)
        push!(w, l)
        if r != l
            push!(pos, s)
            push!(w, r)
        end
    end
    return pos, w
end

"""Value of a profile just before (`left = true`) or just after `s`."""
function _profile_limit(pos::Vector{Float64}, w::Vector{Float64}, s::Float64, left::Bool)
    n = length(pos)
    if left
        i = searchsortedfirst(pos, s)
        i > n && return w[end]
        (pos[i] == s || i == 1) && return w[i]
        return w[i - 1] + (w[i] - w[i - 1]) * (s - pos[i - 1]) / (pos[i] - pos[i - 1])
    end
    i = searchsortedlast(pos, s)
    i < 1 && return w[1]
    (pos[i] == s || i == n) && return w[i]
    return w[i] + (w[i + 1] - w[i]) * (s - pos[i]) / (pos[i + 1] - pos[i])
end

"""Plan area of `outer` minus the parts of `holes` inside it (m²)."""
function _net_area(outer::Vector{NTuple{2,Float64}}, holes::Vector{Vector{NTuple{2,Float64}}})
    lo = minimum(p[1] for p in outer)
    hi = maximum(p[1] for p in outer)
    _, _, inside, _, _ = _slice_profile(outer, holes, lo, hi)
    return inside
end

"""`pts` in the frame of edge `a → b`: `s` from 0 at `a` to 1 at `b`, `d` in meters to the left."""
function _to_edge_frame(pts::Vector{NTuple{2,Float64}}, a::NTuple{2,Float64}, b::NTuple{2,Float64})
    ex, ey = b[1] - a[1], b[2] - a[2]
    L = hypot(ex, ey)
    tx, ty = ex / L, ey / L
    return NTuple{2,Float64}[(((p[1] - a[1]) * tx + (p[2] - a[2]) * ty) / L,
                              -(p[1] - a[1]) * ty + (p[2] - a[2]) * tx) for p in pts]
end

"""
Convert parametric (s,d) to absolute (x,y) coordinates in meters.

beam_start and beam_end must be in meters.
"""
function vertices(trib::TributaryPolygon, beam_start::NTuple{2,Float64}, 
                  beam_end::NTuple{2,Float64})::Vector{NTuple{2,Float64}}
    beam_vec = (beam_end[1] - beam_start[1], beam_end[2] - beam_start[2])
    beam_len = hypot(beam_vec...)
    beam_len < 1e-12 && return NTuple{2, Float64}[]
    
    beam_dir = (beam_vec[1] / beam_len, beam_vec[2] / beam_len)
    beam_normal = (-beam_dir[2], beam_dir[1])
    
    return [(beam_start[1] + s * beam_len * beam_dir[1] + d * beam_normal[1],
             beam_start[2] + s * beam_len * beam_dir[2] + d * beam_normal[2])
            for (s, d) in zip(trib.s, trib.d)]
end

"""Convert Meshes.Point to (x,y) tuple in meters."""
function _to_2d(p::Point)
    c = coords(p)
    (Float64(ustrip(u"m", c.x)), Float64(ustrip(u"m", c.y)))
end

"""Convert absolute vertices to parametric (s,d) relative to beam."""
function _to_parametric(abs_verts::Vector{NTuple{2,Float64}}, beam_start::NTuple{2,Float64},
                        beam_end::NTuple{2,Float64})::Tuple{Vector{Float64}, Vector{Float64}}
    isempty(abs_verts) && return (Float64[], Float64[])
    
    beam_vec = (beam_end[1] - beam_start[1], beam_end[2] - beam_start[2])
    beam_len = hypot(beam_vec...)
    beam_len < 1e-12 && return (zeros(length(abs_verts)), zeros(length(abs_verts)))
    
    beam_dir = (beam_vec[1] / beam_len, beam_vec[2] / beam_len)
    beam_normal = (-beam_dir[2], beam_dir[1])
    
    s_vals, d_vals = Float64[], Float64[]
    for v in abs_verts
        rel = (v[1] - beam_start[1], v[2] - beam_start[2])
        push!(s_vals, (rel[1] * beam_dir[1] + rel[2] * beam_dir[2]) / beam_len)
        push!(d_vals, rel[1] * beam_normal[1] + rel[2] * beam_normal[2])
    end
    
    _rotate_to_beam_first(s_vals, d_vals)
end

"""Rotate (s,d) arrays so beam edge vertices come first."""
function _rotate_to_beam_first(s::Vector{Float64}, d::Vector{Float64})
    n = length(s)
    n < 2 && return (s, d)
    
    tol = 1e-6
    best_idx, best_score = 1, Inf
    for i in 1:n
        j = mod1(i + 1, n)
        score = abs(d[i]) + abs(d[j])
        if score < best_score && s[i] <= s[j] + tol
            best_score, best_idx = score, i
        end
    end
    
    best_idx == 1 && return (s, d)
    (vcat(s[best_idx:end], s[1:best_idx-1]), vcat(d[best_idx:end], d[1:best_idx-1]))
end

# =============================================================================
# Cell Depth Computation (for EFM l2 calculation)
# =============================================================================

"""
    ray_segment_intersection(origin, direction, seg_start, seg_end) -> Union{Float64, Nothing}

Compute the distance t along a ray where it intersects a line segment.

Ray: P(t) = origin + t * direction, t ≥ 0
Returns t if intersection exists (positive direction), nothing otherwise.
"""
function ray_segment_intersection(
    origin::NTuple{2,Float64}, 
    direction::NTuple{2,Float64},
    seg_start::NTuple{2,Float64}, 
    seg_end::NTuple{2,Float64}
)::Union{Float64, Nothing}
    # Ray: P = origin + t * direction
    # Segment: Q = seg_start + u * (seg_end - seg_start), u ∈ [0,1]
    # Solve: origin + t * direction = seg_start + u * seg_vec
    
    seg_vec = (seg_end[1] - seg_start[1], seg_end[2] - seg_start[2])
    
    # Cross product of direction and seg_vec (2D: returns scalar)
    denom = direction[1] * seg_vec[2] - direction[2] * seg_vec[1]
    
    # Parallel lines (no intersection or infinite intersections)
    abs(denom) < 1e-12 && return nothing
    
    # Vector from ray origin to segment start
    diff = (seg_start[1] - origin[1], seg_start[2] - origin[2])
    
    # Solve for t and u
    t = (diff[1] * seg_vec[2] - diff[2] * seg_vec[1]) / denom
    u = (diff[1] * direction[2] - diff[2] * direction[1]) / denom
    
    # Check if intersection is valid (t > 0 and u ∈ [0,1])
    if t > 1e-9 && u >= -1e-9 && u <= 1.0 + 1e-9
        return t
    end
    
    return nothing
end

"""
    _point_in_polygon(pt, poly_pts) -> Bool

Check if a point is strictly inside a polygon using ray casting (crossing number).
"""
function _point_in_polygon(pt::NTuple{2,Float64}, poly_pts::Vector{NTuple{2,Float64}})::Bool
    n = length(poly_pts)
    n < 3 && return false
    
    inside = false
    j = n
    for i in 1:n
        xi, yi = poly_pts[i]
        xj, yj = poly_pts[j]
        
        if ((yi > pt[2]) != (yj > pt[2])) &&
           (pt[1] < (xj - xi) * (pt[2] - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

"""
    perpendicular_depth_at_s(cell_pts, edge_idx, s; debug=false) -> Float64

Compute the perpendicular depth from a point on the edge to the cell boundary.

Uses the "offset origin" technique: start slightly inside the polygon to avoid
boundary degeneracies, then find the first intersection (which is where we exit).

# Arguments
- `cell_pts`: Cell polygon vertices as (x,y) tuples in meters
- `edge_idx`: Index of the edge (1-indexed)
- `s`: Parametric position along edge (0 = start, 1 = end)
- `debug`: If true, print debug information

# Returns
Perpendicular distance to the closest cell boundary in the inward direction (meters).
Returns 0 if the perpendicular direction doesn't point into the polygon interior.
"""
function perpendicular_depth_at_s(
    cell_pts::Vector{NTuple{2,Float64}}, 
    edge_idx::Int, 
    s::Float64;
    debug::Bool = false
)::Float64
    # NOTE: Caller must ensure cell_pts is CCW oriented!
    # We don't call _ensure_ccw here because edge_idx is tied to the input ordering.
    n = length(cell_pts)
    n < 3 && return 0.0
    
    # Get edge endpoints
    v1 = cell_pts[edge_idx]
    v2 = cell_pts[mod1(edge_idx + 1, n)]
    
    # Point on edge at parametric position s
    origin_on_edge = (v1[1] + s * (v2[1] - v1[1]), v1[2] + s * (v2[2] - v1[2]))
    
    # Edge direction and inward perpendicular
    edge_vec = (v2[1] - v1[1], v2[2] - v1[2])
    edge_len = hypot(edge_vec...)
    edge_len < 1e-12 && return 0.0
    
    edge_dir = (edge_vec[1] / edge_len, edge_vec[2] / edge_len)
    # Perpendicular pointing "left" (inward for CCW polygon)
    perp_dir = (-edge_dir[2], edge_dir[1])
    
    # Offset origin slightly inside to avoid boundary degeneracies
    # This is the standard robust trick for ray casting from polygon boundaries
    ε = 1e-6 * edge_len
    origin = (origin_on_edge[1] + ε * perp_dir[1], origin_on_edge[2] + ε * perp_dir[2])
    
    # Check if the offset origin is actually inside the polygon
    # If not, perp_dir is pointing outward (concave case)
    if !_point_in_polygon(origin, cell_pts)
        debug && println("  Edge $edge_idx, s=$s: offset origin outside polygon")
        return 0.0
    end
    
    if debug
        println("  Edge $edge_idx, s=$s: origin=$origin_on_edge, perp_dir=$perp_dir")
    end
    
    # Find closest intersection with any OTHER edge
    # We check all edges except the current one; the point-in-polygon check above
    # already ensures we're starting inside, so the first hit is valid
    min_t = Inf
    best_edge = -1
    
    for i in 1:n
        # Only skip the edge we're measuring from
        i == edge_idx && continue
        
        seg_start = cell_pts[i]
        seg_end = cell_pts[mod1(i + 1, n)]
        
        t = ray_segment_intersection(origin, perp_dir, seg_start, seg_end)
        if !isnothing(t) && t > 0 && t < min_t
            min_t = t
            best_edge = i
            debug && println("    Hit edge $i at t=$t")
        end
    end
    
    if debug
        depth = isfinite(min_t) ? (min_t + ε) : 0.0
        println("    -> Best: edge $best_edge, depth=$depth")
    end
    
    # Add back the ε offset to get true distance from edge
    return isfinite(min_t) ? (min_t + ε) : 0.0
end

"""
    compute_cell_depths(cell_pts, edge_idx) -> NamedTuple

Compute cell boundary depths for an edge using vertex projection.

Projects all polygon vertices onto the edge to find critical points where
depth changes, then computes the stiffness-consistent effective width (l2_stiff)
via cubic mean integration.

# Arguments
- `cell_pts`: Cell polygon vertices as (x,y) tuples in meters
- `edge_idx`: Index of the edge (1-indexed)

# Returns
NamedTuple with:
- `s_vals`: Parametric positions where depths were sampled
- `depths`: Perpendicular depths at each sample point (m)
- `l2_stiff`: Stiffness-consistent width via cubic mean (m)
- `max`: Maximum depth (m)
- `s_max`: Parametric position of max depth
"""
function compute_cell_depths(
    cell_pts::Vector{NTuple{2,Float64}}, 
    edge_idx::Int
)
    # NOTE: Caller must ensure cell_pts is CCW oriented!
    n = length(cell_pts)
    empty_result = (s_vals=Float64[], depths=Float64[], l2_stiff=0.0, max=0.0, s_max=0.0)
    n < 3 && return empty_result
    
    # Get edge endpoints and properties
    v1 = cell_pts[edge_idx]
    v2 = cell_pts[mod1(edge_idx + 1, n)]
    
    edge_vec = (v2[1] - v1[1], v2[2] - v1[2])
    edge_len = hypot(edge_vec...)
    edge_len < 1e-12 && return empty_result
    
    edge_dir = (edge_vec[1] / edge_len, edge_vec[2] / edge_len)
    
    # =========================================================================
    # Step 1: Find critical s values by projecting vertices onto edge
    # =========================================================================
    # Use open interval ]ε, 1-ε[ to avoid vertex degeneracies
    ϵ = 1e-3  # Stay 0.1% away from endpoints
    
    critical_s = Float64[ϵ, 1.0 - ϵ]  # Near-endpoints, not exact
    
    for i in 1:n
        i == edge_idx && continue  # Skip edge's own vertices
        i == mod1(edge_idx + 1, n) && continue
        
        vtx = cell_pts[i]
        rel = (vtx[1] - v1[1], vtx[2] - v1[2])
        s_proj = (rel[1] * edge_dir[1] + rel[2] * edge_dir[2]) / edge_len
        
        # Clamp to open interval and only include if interior
        s_proj = clamp(s_proj, ϵ, 1.0 - ϵ)
        if s_proj > ϵ + 1e-9 && s_proj < 1.0 - ϵ - 1e-9
            push!(critical_s, s_proj)
        end
    end
    
    # Remove duplicates and sort
    unique!(sort!(critical_s))
    
    # =========================================================================
    # Step 2: Compute depths at all critical points
    # =========================================================================
    depths = [perpendicular_depth_at_s(cell_pts, edge_idx, s) for s in critical_s]
    
    # =========================================================================
    # Step 3: Find max depth and its position
    # =========================================================================
    max_idx = argmax(depths)
    depth_max = depths[max_idx]
    s_max = critical_s[max_idx]
    
    # =========================================================================
    # Step 4: Compute stiffness-consistent L2 via cubic mean
    # L2_stiff = (∫ L(s)³ ds / ∫ ds)^(1/3)
    # This weights longer reaches more heavily - conservative for deflection
    # Since depth is piecewise-linear, trapezoidal rule on L³ is exact per segment
    # =========================================================================
    s_range = critical_s[end] - critical_s[1]
    cubic_integral = 0.0
    for i in 1:(length(critical_s) - 1)
        Δs = critical_s[i+1] - critical_s[i]
        # Trapezoidal rule on L³
        avg_cubic = (depths[i]^3 + depths[i+1]^3) / 2
        cubic_integral += avg_cubic * Δs
    end
    l2_stiff = s_range > 0 ? (cubic_integral / s_range)^(1/3) : 0.0
    
    return (s_vals=critical_s, depths=depths, l2_stiff=l2_stiff, max=depth_max, s_max=s_max)
end

"""
Create TributaryPolygon with cell depth values computed.

This version computes the cell boundary depths (for EFM l2 calculation).

`profile` is the edge's exact `(positions, widths)` when the partition knows
it (one-way). Otherwise the region `verts` minus `holes` is sliced across the
edge ([`_edge_slice_profile`](@ref)), and with holes `area` is replaced by that
net area.
"""
function _make_tributary_with_depths(edge_idx::Int, verts::Vector{NTuple{2,Float64}}, 
                                     original_pts::Vector{NTuple{2,Float64}}, area::Float64, 
                                     frac::Float64;
                                     profile = nothing,
                                     holes::Vector{Vector{NTuple{2,Float64}}} = Vector{NTuple{2,Float64}}[]
                                     )::TributaryPolygon
    n = length(original_pts)
    a = original_pts[edge_idx]
    b = original_pts[mod1(edge_idx + 1, n)]
    
    # Compute parametric representation
    if isempty(verts)
        s, d = Float64[], Float64[]
    else
        s, d = _to_parametric(verts, a, b)
    end

    if profile !== nothing
        ps, pw = profile
    elseif isempty(verts)
        ps, pw = [0.0, 1.0], [0.0, 0.0]
    else
        frame_holes = [_to_edge_frame(h, a, b) for h in holes]
        ps, pw, net = _edge_slice_profile(_to_edge_frame(verts, a, b), frame_holes)
        isempty(holes) || (area = net * hypot(b[1] - a[1], b[2] - a[2]))
    end
    
    # Compute cell boundary depths (includes l2_stiff via cubic mean)
    cell = compute_cell_depths(original_pts, edge_idx)
    
    TributaryPolygon(edge_idx, s, d, area, frac,
                     cell.s_vals, cell.depths, cell.l2_stiff,
                     cell.max, cell.s_max, ps, pw)
end

"""
    assert_tributary_balance(tribs, pts, net_area; context="")

Throw unless every edge's load profile carries its area and the areas sum to
`net_area`. A shortfall is floor load that would reach no beam.
"""
function assert_tributary_balance(tribs::Vector{TributaryPolygon},
                                  pts::Vector{NTuple{2,Float64}}, net_area::Float64;
                                  context::AbstractString = "")
    n = length(pts)
    tol = 1e-6 * max(net_area, 1.0)
    total = 0.0
    for t in tribs
        i = t.local_edge_idx
        L = hypot(pts[mod1(i + 1, n)][1] - pts[i][1], pts[mod1(i + 1, n)][2] - pts[i][2])
        carried = profile_area(t, L)
        abs(carried - t.area) <= tol || error(
            "Tributary of edge $i$(context) has area $(t.area) m² but its load profile " *
            "carries $carried m². The rest of that floor would reach no beam.")
        total += t.area
    end
    abs(total - net_area) <= tol || error(
        "Tributaries$(context) sum to $total m² on a floor of $net_area m². " *
        "The partition " * (total < net_area ? "misses" : "double counts") *
        " $(abs(total - net_area)) m².")
    return nothing
end

"""Shoelace formula for signed polygon area."""
function _polygon_area(pts::Vector{NTuple{2,Float64}})
    n = length(pts)
    n < 3 && return 0.0
    sum(pts[i][1] * pts[mod1(i+1,n)][2] - pts[mod1(i+1,n)][1] * pts[i][2] for i in 1:n) / 2
end

_is_ccw(pts::Vector{NTuple{2,Float64}}) = _polygon_area(pts) > 0

"""Ensure polygon is CCW oriented."""
_ensure_ccw(pts::Vector{NTuple{2,Float64}}) = _is_ccw(pts) ? pts : reverse(pts)

"""Simplify polygon by dropping collinear vertices."""
function simplify_collinear_polygon(pts::Vector{NTuple{2,Float64}}; tol=1e-12)
    n = length(pts)
    n ≤ 3 && return pts, collect(1:n)
    
    is_collinear(i) = begin
        p_prev, p, p_next = pts[mod1(i-1,n)], pts[i], pts[mod1(i+1,n)]
        abs((p[1]-p_prev[1])*(p_next[2]-p[2]) - (p[2]-p_prev[2])*(p_next[1]-p[1])) ≤ tol
    end
    
    keep = [i for i in 1:n if !is_collinear(i)]
    simp = [pts[i] for i in keep]
    
    if _polygon_area(simp) < 0
        reverse!(simp)
        reverse!(keep)
    end
    
    simp, keep
end
