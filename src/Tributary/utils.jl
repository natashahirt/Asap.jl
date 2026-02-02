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
end

# Constructor with default cell depths (for backwards compatibility)
function TributaryPolygon(local_edge_idx::Int, s::Vector{Float64}, d::Vector{Float64}, 
                          area::Float64, fraction::Float64)
    TributaryPolygon(local_edge_idx, s, d, area, fraction, 
                     Float64[], Float64[], 0.0, 0.0, 0.0)
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
"""
function _make_tributary_with_depths(edge_idx::Int, verts::Vector{NTuple{2,Float64}}, 
                                     original_pts::Vector{NTuple{2,Float64}}, area::Float64, 
                                     frac::Float64)::TributaryPolygon
    n = length(original_pts)
    
    # Compute parametric representation
    if isempty(verts)
        s, d = Float64[], Float64[]
    else
        s, d = _to_parametric(verts, original_pts[edge_idx], original_pts[mod1(edge_idx + 1, n)])
    end
    
    # Compute cell boundary depths (includes l2_stiff via cubic mean)
    cell = compute_cell_depths(original_pts, edge_idx)
    
    TributaryPolygon(edge_idx, s, d, area, frac,
                     cell.s_vals, cell.depths, cell.l2_stiff,
                     cell.max, cell.s_max)
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
