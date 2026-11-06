# ============================================================================
# Structured Mesh Generation (adapted from FinEtools.jl MeshTriangleModule)
# ============================================================================

"""
    _t3block(Lx, Ly, nx, ny) -> (points, connectivity)

Generate a structured triangular mesh of a rectangle [0,Lx] × [0,Ly].

Returns `(points, conn)` where:
- `points::Vector{Tuple{Float64,Float64}}` — node coordinates
- `conn::Vector{NTuple{3,Int}}` — triangle connectivity (1-indexed)

Each rectangular cell is split into 2 triangles (lower-left diagonal).
Produces `(nx+1)*(ny+1)` nodes and `2*nx*ny` triangles.

Adapted from FinEtools.jl `T3block` by Petr Krysl (MIT License).
"""
function _t3block(Lx::Float64, Ly::Float64, nx::Int, ny::Int)
    xs = range(0.0, Lx, length=nx+1)
    ys = range(0.0, Ly, length=ny+1)
    return _t3blockx(collect(xs), collect(ys))
end

"""
    _t3blockx(xs, ys) -> (points, connectivity)

Generate a structured triangular mesh from explicit coordinate arrays.
Supports non-uniform (graded) spacing.

**Alternating diagonals** (checkerboard pattern): each quad cell `(i,j)` is
split along the `/` diagonal when `(i+j)` is even and along the `\\` diagonal
when `(i+j)` is odd. This eliminates the directional stiffness bias that a
uniform diagonal orientation would introduce in FE analysis.

Returns `(points, conn)` where:
- `points::Vector{Tuple{Float64,Float64}}` — node coordinates
- `conn::Vector{NTuple{3,Int}}` — triangle connectivity (1-indexed)

Adapted from FinEtools.jl `T3blockx` by Petr Krysl (MIT License), with
alternating-diagonal enhancement for structural analysis.
"""
function _t3blockx(xs::Vector{Float64}, ys::Vector{Float64})
    nL = length(xs) - 1
    nW = length(ys) - 1
    nnodes = (nL + 1) * (nW + 1)
    ncells = 2 * nL * nW

    points = Vector{Tuple{Float64, Float64}}(undef, nnodes)
    conn = Vector{NTuple{3, Int}}(undef, ncells)

    # Generate node coordinates
    f = 1
    for j in 1:(nW+1)
        for i in 1:(nL+1)
            points[f] = (xs[i], ys[j])
            f += 1
        end
    end

    # Generate triangle connectivity with alternating diagonals
    # Node numbering for cell (i, j):
    #   nw = f+(nL+1)   ne = f+(nL+1)+1
    #   sw = f           se = f+1
    gc = 1
    for i in 1:nL
        for j in 1:nW
            sw = (j - 1) * (nL + 1) + i
            se = sw + 1
            nw = sw + (nL + 1)
            ne = nw + 1
            
            if iseven(i + j)
                # "/" diagonal: sw → ne
                conn[gc]     = (sw, se, ne)
                conn[gc + 1] = (sw, ne, nw)
            else
                # "\" diagonal: se → nw
                conn[gc]     = (sw, se, nw)
                conn[gc + 1] = (se, ne, nw)
            end
            gc += 2
        end
    end

    return points, conn
end

# ============================================================================
# Graded Spacing for Mesh Refinement
# ============================================================================

"""
    _graded_spacing(x0, x1, targets, h_far, h_near, r_transition) -> Vector{Float64}

Build a non-uniform coordinate array that clusters points near `targets`.
Spacing is **symmetric** around each target — the algorithm walks from both ends
of every inter-target segment toward the midpoint and merges the results.

- `x0, x1`: domain bounds
- `targets`: locations requiring refinement (e.g., column coordinates)
- `h_far`: far-field element size
- `h_near`: near-target element size (must be ≤ h_far)
- `r_transition`: distance over which spacing grades from h_near to h_far
"""
function _graded_spacing(x0::Float64, x1::Float64, targets::Vector{Float64},
                         h_far::Float64, h_near::Float64, r_transition::Float64)
    @assert h_near <= h_far "h_near ($h_near) must be ≤ h_far ($h_far)"
    @assert x0 < x1 "x0 must be < x1"
    @assert r_transition > 0 "r_transition must be > 0"
    
    # Protected target points (must survive dedup)
    target_set = Set{Float64}([x0, x1])
    for t in targets
        if x0 ≤ t ≤ x1
            push!(target_set, t)
        end
    end
    
    # Sorted list of mandatory "waypoints" (endpoints + interior targets)
    waypoints = sort!(collect(target_set))
    
    # Collect all points
    pts = Set{Float64}()
    for wp in waypoints
        push!(pts, wp)
    end
    
    # Helper: local element size at position x
    _h_at(x) = begin
        d_min = Inf
        for t in waypoints
            d_min = min(d_min, abs(x - t))
        end
        ratio = clamp(d_min / r_transition, 0.0, 1.0)
        h_near + ratio * (h_far - h_near)
    end
    
    # For each segment [a, b], walk from BOTH ends toward the midpoint.
    # This ensures spacing is symmetric around each waypoint.
    for seg in 1:length(waypoints)-1
        a, b = waypoints[seg], waypoints[seg+1]
        mid = (a + b) / 2.0
        
        # Forward walk: a → mid
        x = a
        while x < mid - h_near * 0.05
            h_local = _h_at(x)
            x_next = min(x + h_local, mid)
            push!(pts, x_next)
            x = x_next
        end
        
        # Backward walk: b → mid
        x = b
        while x > mid + h_near * 0.05
            h_local = _h_at(x)
            x_next = max(x - h_local, mid)
            push!(pts, x_next)
            x = x_next
        end
    end
    
    result = sort!(collect(pts))
    
    # Dedup close points, but never remove protected targets
    min_spacing = h_near * 0.3
    filtered = Float64[result[1]]
    for i in 2:length(result)
        gap = result[i] - filtered[end]
        is_protected = result[i] in target_set
        if gap > min_spacing || i == length(result) || is_protected
            push!(filtered, result[i])
        end
    end
    if filtered[end] != x1
        push!(filtered, x1)
    end
    
    return filtered
end

"""
    _warn_mesh_density(h, Lx, Ly; reference_shortest_side=nothing)

Warn if the effective mesh density seems unreasonably coarse or fine.
Thresholds are expressed in terms of effective `n` (divisions along a reference shortest side).

When `reference_shortest_side` is set (meters), it is used instead of `min(Lx, Ly)` from the
panel bounding box. Callers meshing a **multi-cell** slab as one polygon should pass the
typical cell shortest span so `n` reflects per-cell density (matching how `target_edge_length`
is often chosen), not divisions across the full outline.
"""
function _warn_mesh_density(h::Float64, Lx::Float64, Ly::Float64;
                          reference_shortest_side::Union{Float64, Nothing} = nothing)
    L_min = something(reference_shortest_side, min(Lx, Ly))
    effective_n = L_min / h
    if effective_n < 4
        @warn "Shell mesh is very coarse (effective n ≈ $(round(effective_n, digits=1)) " *
              "per reference shortest side). Consider reducing target_edge_length."
    elseif effective_n > 100
        @warn "Shell mesh is very fine (effective n ≈ $(round(effective_n, digits=1)) " *
              "per reference shortest side). Consider increasing target_edge_length for faster analysis."
    end
end

"""
    _uniform_spacing(x0, x1, h) -> Vector{Float64}

Build a uniform coordinate array with target spacing `h`.
"""
function _uniform_spacing(x0::Float64, x1::Float64, h::Float64)
    n = max(1, ceil(Int, (x1 - x0) / h))
    return collect(range(x0, x1, length=n+1))
end

"""
    _is_rectangular(boundary_pts; tol=1e-6) -> Bool

Check if a polygon is a (possibly rotated) axis-aligned rectangle.
"""
function _is_rectangular(boundary_pts::Vector{Tuple{Float64, Float64}}; tol::Float64=1e-6)
    length(boundary_pts) != 4 && return false
    
    # Check that opposite edges are parallel and equal length
    p1, p2, p3, p4 = boundary_pts
    
    # Edge vectors
    e1 = (p2[1] - p1[1], p2[2] - p1[2])
    e2 = (p3[1] - p2[1], p3[2] - p2[2])
    e3 = (p4[1] - p3[1], p4[2] - p3[2])
    e4 = (p1[1] - p4[1], p1[2] - p4[2])
    
    # Check axis-aligned: each edge should be horizontal or vertical
    h1 = abs(e1[2]) < tol  # horizontal
    v1 = abs(e1[1]) < tol  # vertical
    h2 = abs(e2[2]) < tol
    v2 = abs(e2[1]) < tol
    h3 = abs(e3[2]) < tol
    v3 = abs(e3[1]) < tol
    h4 = abs(e4[2]) < tol
    v4 = abs(e4[1]) < tol
    
    return (h1 || v1) && (h2 || v2) && (h3 || v3) && (h4 || v4)
end

