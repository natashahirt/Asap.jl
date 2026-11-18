# =============================================================================
# Polygon Span Calculations
# =============================================================================
#
# A span is the longest straight segment that stays inside the slab region:
# the outer boundary minus any voids. A directional span fixes the segment
# direction; the isotropic span allows any direction.
#
# A line that only touches the boundary at a corner is measured as the limit of
# the parallel lines beside it. A span therefore never runs along a zero-width
# strip pinched between two notches, and a corner on a supported boundary ends
# a span unless the strip beside it runs clear.

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
- `isotropic::T`: Longest internal segment in any direction

## Notes
- For one-way systems, use `primary` for sizing
- For two-way systems, use [`longer_axis_span`](@ref): the longer of the two
  reinforcement-axis spans. `isotropic` is the longest segment in any
  direction, which is the diagonal of a rectangle, not a reinforcement span
- `axis` is the primary reinforcement direction. With no given axis it follows
  the slab's orthogonal edge frame, and falls back to global x when the edges
  do not share one frame
"""
struct SpanInfo{T}
    primary::T
    secondary::T
    axis::NTuple{2, Float64}
    isotropic::T
end

"""
    SpanInfo(verts; axis=nothing, holes=nothing)

Compute SpanInfo from polygon vertices.

## Arguments
- `verts::Vector{<:Point}`: Polygon vertices as Meshes.Point objects
- `axis::Union{Nothing, NTuple{2,Float64}}`: Optional primary span direction.
  If `nothing`, the short axis of the slab's orthogonal edge frame. A bay
  whose edges do not share one frame (a skewed parallelogram) keeps global x.
- `holes`: Optional voids, each a vertex loop (`Point`s or `(x, y)` in meters).

## Returns
- `SpanInfo{Float64}`: Span info with values in meters (unitless Float64)
"""
function SpanInfo(verts::Vector{<:Point};
                  axis::Union{Nothing, NTuple{2,Float64}}=nothing,
                  holes=nothing)
    outer, hs = _span_region(verts, holes)
    iso = _span_isotropic(outer, hs)

    if isnothing(axis)
        # Short axis of the orthogonal frame the edges agree on. Global x
        # when they do not, which is the historical span pair.
        frame = _orthogonal_frame_axis(outer)
        ax = something(frame, (1.0, 0.0))
        perp = (-ax[2], ax[1])
        along = _span_directional(outer, hs, ax)
        across = _span_directional(outer, hs, perp)
        if along <= across
            primary, secondary = along, across
            axis = ax
        else
            primary, secondary = across, along
            axis = perp
        end
    else
        ax = (Float64(axis[1]), Float64(axis[2]))
        alen = hypot(ax...)
        alen < 1e-12 && error("axis must be non-zero")
        ax = (ax[1] / alen, ax[2] / alen)
        perp = (-ax[2], ax[1])

        primary = _span_directional(outer, hs, ax)
        secondary = _span_directional(outer, hs, perp)
        axis = ax
    end

    SpanInfo{Float64}(primary, secondary, axis, iso)
end

"""
    governing_spans(span_infos::Vector{SpanInfo{T}}) -> SpanInfo{T}

Governing spans of several cells that share one reinforcement frame.

Each cell's spans are taken along the first cell's axis and the perpendicular,
so a cell whose primary axis is turned 90° still contributes its span to the
right direction. Cells whose axes are neither parallel nor perpendicular do
not share a frame, and this throws rather than pairing unrelated spans.

`primary` stays the span along the first cell's axis. It is not re-sorted to
be the shorter one.
"""
function governing_spans(span_infos::Vector{SpanInfo{T}}) where T
    isempty(span_infos) && error("Cannot compute governing spans from empty vector")
    length(span_infos) == 1 && return span_infos[1]

    ref = span_infos[1].axis
    along = span_infos[1].primary
    across = span_infos[1].secondary
    isotropic = span_infos[1].isotropic
    for si in @view span_infos[2:end]
        d = abs(si.axis[1] * ref[1] + si.axis[2] * ref[2])
        if d >= 1 - 1e-6
            along = max(along, si.primary)
            across = max(across, si.secondary)
        elseif d <= 1e-6
            along = max(along, si.secondary)
            across = max(across, si.primary)
        else
            throw(ArgumentError(
                "cells in one slab do not share a reinforcement axis: " *
                "$(ref) and $(si.axis)"))
        end
        isotropic = max(isotropic, si.isotropic)
    end
    return SpanInfo{T}(along, across, ref, isotropic)
end

# Convenience accessors
short_span(si::SpanInfo) = si.primary
long_span(si::SpanInfo) = si.secondary
two_way_span(si::SpanInfo) = si.isotropic

"""
    longer_axis_span(si) -> span

Longer of the two reinforcement-axis spans. Two-way minimum thickness uses
this where a code table says "the longer span". It is not [`two_way_span`](@ref),
which is the longest segment in any direction, and it is not the clear
distance between column faces.
"""
longer_axis_span(si::SpanInfo) = max(si.primary, si.secondary)

"""
Unit axis of the orthogonal frame the boundary edges agree on, or `nothing`.

Edge directions fold into a right angle, because a reinforcement pair repeats
every 90°, and each edge votes with its length. A frame within half a degree
of the global axes snaps to `(1, 0)`, so an axis-aligned bay keeps that exact
tuple. Edges that do not agree (a skewed parallelogram, a fan of curve
facets) return `nothing` and the caller keeps the global axes.
"""
function _orthogonal_frame_axis(outer)::Union{Nothing, NTuple{2, Float64}}
    length(outer) < 3 && return nothing
    s = 0.0
    c = 0.0
    w = 0.0
    n = length(outer)
    @inbounds for i in 1:n
        p, q = outer[i], outer[mod1(i + 1, n)]
        dx, dy = q[1] - p[1], q[2] - p[2]
        L = hypot(dx, dy)
        L <= 1e-9 && continue
        α = mod(atan(dy, dx), π / 2)
        s += L * sin(4α)
        c += L * cos(4α)
        w += L
    end
    w <= 1e-9 && return nothing
    hypot(s, c) / w < 0.95 && return nothing
    θ = mod(atan(s, c), 2π) / 4
    min(θ, π / 2 - θ) < deg2rad(0.5) && return (1.0, 0.0)
    return (cos(θ), sin(θ))
end

# =============================================================================
# Public Span Function
# =============================================================================

"""
    get_polygon_span(vertices; axis=nothing, holes=nothing)

Longest straight segment inside the polygon, not crossing any void.

## Arguments
- `vertices::Vector{<:Point}`: Outer boundary as Meshes.Point objects
- `axis::Union{Nothing, AbstractVector{<:Real}}`: Direction vector [vx, vy].
  - If `nothing` (or zero), returns the **isotropic span**: longest segment in any direction
  - If provided, returns the **directional span**: longest segment parallel to `axis`
- `holes`: Optional voids, each a vertex loop (`Point`s or `(x, y)` in meters).
  Voids may overlap the outer boundary.

## Returns
- `Float64`: Maximum span in meters

## Method
Every boundary edge is straight, so as a line slides parallel to itself the
length of each inside piece changes linearly between vertices. The longest
directional piece therefore lies on a line through a vertex, taken as the limit
from either side. Rotating a line about a vertex makes each piece length convex
in the angle, so the longest isotropic piece lies on a line through two vertices.

## Examples
```julia
span = get_polygon_span(vertices)                       # isotropic
span_x = get_polygon_span(vertices; axis=[1.0, 0.0])    # along x
span_y = get_polygon_span(vertices; axis=[0.0, 1.0], holes=[opening])
```
"""
function get_polygon_span(
    vertices::Vector{<:Point};
    axis::Union{Nothing, AbstractVector{<:Real}} = nothing,
    holes = nothing,
)
    outer, hs = _span_region(vertices, holes)
    if isnothing(axis) || hypot(axis[1], axis[2]) < 1e-12
        return _span_isotropic(outer, hs)
    end
    vlen = hypot(axis[1], axis[2])
    return _span_directional(outer, hs, (Float64(axis[1]) / vlen, Float64(axis[2]) / vlen))
end

# =============================================================================
# Region Setup
# =============================================================================

_span_pt(p::Point) = _to_2d(p)
_span_pt(p) = (Float64(p[1]), Float64(p[2]))

"""Vertex loop in meters, without a repeated closing vertex."""
function _span_loop(vertices)
    loop = NTuple{2, Float64}[_span_pt(v) for v in vertices]
    length(loop) > 1 && loop[end] == loop[1] && pop!(loop)
    return loop
end

"""Outer loop and void loops (voids with fewer than 3 vertices are dropped)."""
function _span_region(vertices, holes)
    outer = _span_loop(vertices)
    hs = Vector{NTuple{2, Float64}}[]
    if holes !== nothing
        for h in holes
            loop = _span_loop(h)
            length(loop) >= 3 && push!(hs, loop)
        end
    end
    return outer, hs
end

"""Distance below which a vertex counts as on a line: 1e-9 of the region size."""
function _span_tol(outer)
    xs = (p[1] for p in outer)
    ys = (p[2] for p in outer)
    return 1e-9 * max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys), 1.0)
end

function _max_vertex_distance(pts)
    best = 0.0
    @inbounds for i in 1:length(pts)-1, j in i+1:length(pts)
        best = max(best, hypot(pts[j][1] - pts[i][1], pts[j][2] - pts[i][2]))
    end
    return best
end

"""True if every turn of the loop goes the same way (collinear vertices allowed)."""
function _is_convex(loop, tol)
    m = length(loop)
    m < 4 && return true
    turn = 0
    @inbounds for i in 1:m
        a, b, c = loop[i], loop[mod1(i + 1, m)], loop[mod1(i + 2, m)]
        e1 = (b[1] - a[1], b[2] - a[2])
        e2 = (c[1] - b[1], c[2] - b[2])
        cr = e1[1] * e2[2] - e1[2] * e2[1]
        abs(cr) <= tol * (hypot(e1...) + hypot(e2...)) && continue
        s = cr > 0 ? 1 : -1
        turn == 0 && (turn = s)
        s != turn && return false
    end
    return true
end

# =============================================================================
# Pieces of One Line Inside the Region
# =============================================================================

"""
Which side of the line a vertex lying on it is treated as being on.

Vertices before `split` (by line parameter) go to side `first`; the rest go to
`-first`. A parallel shift puts every on-line vertex on the same side
(`split = Inf`); a small rotation puts vertices on either side of the pivot on
opposite sides.
"""
struct _LineSide
    split::Float64
    first::Int
end
(s::_LineSide)(t::Float64) = t < s.split ? s.first : -s.first

"""Parameter along the line and signed offset from it for point `v`."""
@inline function _line_coords(v, p, d, n)
    dx, dy = v[1] - p[1], v[2] - p[2]
    return dx * d[1] + dy * d[2], dx * n[1] + dy * n[2]
end

"""
Sorted parameters where `loop` crosses the line `p + t·d` (unit normal `n`).

Consecutive pairs bound the inside of the loop (even-odd rule).
"""
function _loop_crossings(loop, p, d, n, tol, side::_LineSide)
    ts = Float64[]
    m = length(loop)
    @inbounds for i in 1:m
        ta, ha = _line_coords(loop[i], p, d, n)
        tb, hb = _line_coords(loop[i == m ? 1 : i + 1], p, d, n)
        on_a, on_b = abs(ha) <= tol, abs(hb) <= tol
        sa = on_a ? side(ta) : (ha > 0 ? 1 : -1)
        sb = on_b ? side(tb) : (hb > 0 ? 1 : -1)
        sa == sb && continue
        t = on_a && on_b ? side.split :
            on_a ? ta :
            on_b ? tb :
            ta + (tb - ta) * ha / (ha - hb)
        push!(ts, t)
    end
    return sort!(ts)
end

"""Longest piece of the line that is inside `outer` and outside every void."""
function _longest_piece(outer, holes, p, d, n, tol, side::_LineSide)
    ts = _loop_crossings(outer, p, d, n, tol, side)
    length(ts) < 2 && return 0.0

    voids = NTuple{2, Float64}[]
    for h in holes
        hts = _loop_crossings(h, p, d, n, tol, side)
        for k in 1:2:length(hts)-1
            push!(voids, (hts[k], hts[k + 1]))
        end
    end
    sort!(voids)

    best = 0.0
    for k in 1:2:length(ts)-1
        lo, hi = ts[k], ts[k + 1]
        cur = lo
        for (a, b) in voids
            b <= cur && continue
            a >= hi && break
            a > cur && (best = max(best, a - cur))
            cur = max(cur, b)
        end
        cur < hi && (best = max(best, hi - cur))
    end
    return best
end

# =============================================================================
# Directional and Isotropic Spans
# =============================================================================

"""Longest segment parallel to unit vector `d` inside the region."""
function _span_directional(outer, holes, d::NTuple{2, Float64})
    length(outer) < 3 && return 0.0
    n = (-d[2], d[1])
    tol = _span_tol(outer)

    offsets = Float64[]
    for loop in (outer, holes...), v in loop
        push!(offsets, v[1] * n[1] + v[2] * n[2])
    end
    sort!(offsets)

    best = 0.0
    last = -Inf
    for h in offsets
        h - last <= tol && continue
        last = h
        p = (h * n[1], h * n[2])
        for s in (1, -1)
            best = max(best, _longest_piece(outer, holes, p, d, n, tol, _LineSide(Inf, s)))
        end
    end
    return best
end

"""Longest segment in any direction inside the region."""
function _span_isotropic(outer, holes)
    pts = isempty(holes) ? outer : vcat(outer, holes...)
    length(outer) < 3 && return _max_vertex_distance(pts)
    tol = _span_tol(outer)
    isempty(holes) && _is_convex(outer, tol) && return _max_vertex_distance(outer)

    best = 0.0
    on_t = Float64[]
    m = length(pts)
    for i in 1:m-1, j in i+1:m
        p, q = pts[i], pts[j]
        len = hypot(q[1] - p[1], q[2] - p[2])
        len <= tol && continue
        d = ((q[1] - p[1]) / len, (q[2] - p[2]) / len)
        n = (-d[2], d[1])

        empty!(on_t)
        for v in pts
            t, h = _line_coords(v, p, d, n)
            abs(h) <= tol && push!(on_t, t)
        end
        sort!(on_t)

        for s in (1, -1)
            best = max(best, _longest_piece(outer, holes, p, d, n, tol, _LineSide(Inf, s)))
            for k in 1:length(on_t)-1
                on_t[k + 1] - on_t[k] <= tol && continue
                split = (on_t[k] + on_t[k + 1]) / 2
                best = max(best, _longest_piece(outer, holes, p, d, n, tol, _LineSide(split, s)))
            end
        end
    end
    return best
end
