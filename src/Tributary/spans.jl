# =============================================================================
# Polygon Span Calculations
# =============================================================================

# =============================================================================
# SpanInfo Type
# =============================================================================

"""
    SpanInfo{T}

Encapsulates span information for a slab bay or cell.

## Fields
- `primary::T`: Governing span along the structural (short) direction
- `secondary::T`: Span perpendicular to primary (long direction)
- `axis::NTuple{2, Float64}`: Unit vector of primary span direction
- `isotropic::T`: Two-way span (max vertex-to-vertex distance)

## Notes
- For one-way systems, use `primary` for sizing
- For two-way systems, use `isotropic` or `max(primary, secondary)`
- `axis` indicates the direction loads span (perpendicular to supporting beams)
"""
struct SpanInfo{T}
    primary::T
    secondary::T
    axis::NTuple{2, Float64}
    isotropic::T
end

"""
    SpanInfo(verts; axis=nothing)

Compute SpanInfo from polygon vertices.

## Arguments
- `verts::Vector{<:Point}`: Polygon vertices as Meshes.Point objects
- `axis::Union{Nothing, NTuple{2,Float64}}`: Optional primary span direction.
  If `nothing`, auto-detects short axis from bounding box.

## Returns
- `SpanInfo{Float64}`: Span info with values in meters (unitless Float64)
"""
function SpanInfo(verts::Vector{<:Point}; axis::Union{Nothing, NTuple{2,Float64}}=nothing)
    iso = get_polygon_span(verts)  # two-way span
    
    if isnothing(axis)
        # Auto-detect: primary = short axis
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        if span_x <= span_y
            primary, secondary = span_x, span_y
            axis = (1.0, 0.0)
        else
            primary, secondary = span_y, span_x
            axis = (0.0, 1.0)
        end
    else
        # Custom axis provided
        ax = (Float64(axis[1]), Float64(axis[2]))
        alen = hypot(ax...)
        alen < 1e-12 && error("axis must be non-zero")
        ax = (ax[1] / alen, ax[2] / alen)
        perp = (-ax[2], ax[1])
        
        primary = get_polygon_span(verts; axis=collect(ax))
        secondary = get_polygon_span(verts; axis=collect(perp))
        axis = ax
    end
    
    SpanInfo{Float64}(primary, secondary, axis, iso)
end

"""
    governing_spans(span_infos::Vector{SpanInfo{T}}) -> SpanInfo{T}

Compute governing SpanInfo from multiple cells (e.g., for a multi-cell slab).
Takes the maximum of each span component across all inputs.

## Notes
- Uses the axis from the first SpanInfo (assumes consistent orientation)
- For mixed orientations, consider computing spans per-axis separately
"""
function governing_spans(span_infos::Vector{SpanInfo{T}}) where T
    isempty(span_infos) && error("Cannot compute governing spans from empty vector")
    length(span_infos) == 1 && return span_infos[1]
    
    primary = maximum(si.primary for si in span_infos)
    secondary = maximum(si.secondary for si in span_infos)
    isotropic = maximum(si.isotropic for si in span_infos)
    axis = span_infos[1].axis  # Use first cell's axis as reference
    
    SpanInfo{T}(primary, secondary, axis, isotropic)
end

# Convenience accessors
short_span(si::SpanInfo) = si.primary
long_span(si::SpanInfo) = si.secondary
two_way_span(si::SpanInfo) = si.isotropic

# =============================================================================
# Low-Level Span Functions
# =============================================================================

"""
    get_polygon_span(vertices; axis=nothing)

Compute the maximum internal span of a polygon.

## Arguments
- `vertices::Vector{<:Point}`: Polygon vertices as Meshes.Point objects
- `axis::Union{Nothing, AbstractVector{<:Real}}`: Direction vector [vx, vy].
  - If `nothing`, returns **two-way span**: max vertex-to-vertex distance (polygon diameter)
  - If provided, returns **one-way span**: max internal chord length along axis direction

## Returns
- `Float64`: Maximum span in meters

## Examples
```julia
# Two-way span (max distance between any two vertices)
span = get_polygon_span(vertices)

# One-way span along x-axis
span_x = get_polygon_span(vertices; axis=[1.0, 0.0])

# One-way span along y-axis
span_y = get_polygon_span(vertices; axis=[0.0, 1.0])
```
"""
function get_polygon_span(
    vertices::Vector{<:Point};
    axis::Union{Nothing, AbstractVector{<:Real}} = nothing
)
    if isnothing(axis) || hypot(axis[1], axis[2]) < 1e-12
        return _get_span_two_way(vertices)
    else
        return _get_span_one_way(vertices, axis)
    end
end

"""
Two-way span: maximum distance between any two vertices.
This is the polygon's diameter (for convex polygons, equal to the farthest pair).
"""
function _get_span_two_way(vertices::Vector{<:Point})
    n = length(vertices)
    n < 2 && return 0.0
    
    pts = [_to_2d(v) for v in vertices]
    
    max_dist = 0.0
    for i in 1:(n-1)
        for j in (i+1):n
            d = hypot(pts[i][1] - pts[j][1], pts[i][2] - pts[j][2])
            max_dist = max(max_dist, d)
        end
    end
    
    return max_dist
end

"""
One-way span: maximum internal chord length along specified axis direction.

For each vertex, casts rays in both directions along the axis and measures
the total internal distance through that point.
"""
function _get_span_one_way(vertices::Vector{<:Point}, axis::AbstractVector{<:Real})
    n = length(vertices)
    n < 3 && return 0.0
    
    pts = [_to_2d(v) for v in vertices]
    pts = _ensure_ccw(pts)
    
    # Normalize axis direction
    vx, vy = Float64(axis[1]), Float64(axis[2])
    vlen = hypot(vx, vy)
    vlen < 1e-12 && return 0.0
    dir = (vx / vlen, vy / vlen)
    neg_dir = (-dir[1], -dir[2])
    
    max_span = 0.0
    
    # Method 1: From each vertex, find chord length through it
    for i in 1:n
        p = pts[i]
        dist_pos = _ray_polygon_distance(p, dir, pts)
        dist_neg = _ray_polygon_distance(p, neg_dir, pts)
        # Chord through vertex = dist in both directions
        chord = dist_pos + dist_neg
        max_span = max(max_span, chord)
    end
    
    # Method 2: Also check edges parallel to axis (their lengths)
    for i in 1:n
        a = pts[i]
        b = pts[mod1(i + 1, n)]
        edge_vec = (b[1] - a[1], b[2] - a[2])
        edge_len = hypot(edge_vec...)
        edge_len < 1e-12 && continue
        
        # Check if edge is approximately parallel to axis
        edge_dir = (edge_vec[1] / edge_len, edge_vec[2] / edge_len)
        dot = abs(edge_dir[1] * dir[1] + edge_dir[2] * dir[2])
        if dot > 0.99  # nearly parallel
            max_span = max(max_span, edge_len)
        end
    end
    
    return max_span
end

"""
Find distance from point p along direction dir to the polygon boundary.
Returns the distance to the first boundary intersection, or 0.0 if none found.
"""
function _ray_polygon_distance(p::NTuple{2,Float64}, dir::NTuple{2,Float64}, 
                                pts::Vector{NTuple{2,Float64}})
    n = length(pts)
    min_t = Inf
    
    for i in 1:n
        a = pts[i]
        b = pts[mod1(i + 1, n)]
        t = _ray_segment_intersect(p, dir, a, b)
        if t > 1e-9 && t < min_t
            min_t = t
        end
    end
    
    return isfinite(min_t) ? min_t : 0.0
end

"""
Compute parameter t where ray (p + t*dir) intersects line segment [a, b].
Returns Inf if no valid intersection.
"""
function _ray_segment_intersect(p::NTuple{2,Float64}, dir::NTuple{2,Float64},
                                 a::NTuple{2,Float64}, b::NTuple{2,Float64})
    # Ray: p + t * dir
    # Segment: a + s * (b - a),  s ∈ [0, 1]
    dx, dy = dir
    ex, ey = b[1] - a[1], b[2] - a[2]
    fx, fy = a[1] - p[1], a[2] - p[2]
    
    denom = dx * ey - dy * ex
    abs(denom) < 1e-12 && return Inf  # parallel
    
    t = (fx * ey - fy * ex) / denom
    s = (fx * dy - fy * dx) / denom
    
    # Valid if on segment and in positive ray direction
    (s >= -1e-9 && s <= 1.0 + 1e-9 && t > 1e-9) ? t : Inf
end
