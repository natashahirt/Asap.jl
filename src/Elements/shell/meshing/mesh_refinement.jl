
# ============================================================================
# Internal helpers
# ============================================================================

"""
Create a Node with the specified fixity.
"""
function _create_node_with_fixity(pos, fixity::Symbol, id::Symbol)
    return Node(pos, fixity, id)
end

function _create_node_with_fixity(pos, fixity::Vector{Bool}, id::Symbol)
    @assert length(fixity) == 6 "Custom fixity must be a 6-element Bool vector"
    node = Node(pos, :free, id)
    node.dof .= fixity
    return node
end

"""
Check if a point is on any boundary edge.
"""
function _is_on_boundary_edge(pt::Tuple{Float64, Float64}, boundary_pts::Vector{Tuple{Float64, Float64}}, tol::Float64)
    nc = length(boundary_pts)
    for i in 1:nc
        p1 = boundary_pts[i]
        p2 = boundary_pts[mod1(i+1, nc)]
        if _point_on_segment(pt, p1, p2, tol)
            return true
        end
    end
    return false
end

"""
Check if a point lies on a line segment.
"""
function _point_on_segment(pt::Tuple{Float64, Float64}, p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, tol::Float64)
    # Check if pt is collinear with p1-p2 and within the segment bounds
    dx, dy = p2[1] - p1[1], p2[2] - p1[2]
    len = hypot(dx, dy)
    if len < tol
        return hypot(pt[1] - p1[1], pt[2] - p1[2]) < tol
    end
    
    # Project pt onto line
    t = ((pt[1] - p1[1]) * dx + (pt[2] - p1[2]) * dy) / (len * len)
    
    # Check if within segment
    if t < -tol/len || t > 1 + tol/len
        return false
    end
    
    # Check distance from line
    proj_x = p1[1] + t * dx
    proj_y = p1[2] + t * dy
    dist = hypot(pt[1] - proj_x, pt[2] - proj_y)
    
    return dist < tol
end

"""
Create shells from triangulation and node map.
"""
function _create_shells_from_tri(
    tri::DT.Triangulation,
    node_map::Dict{Int, Node},
    section::ShellSection,
    id::Symbol
)
    shells = ShellTri3[]
    thickness_q = section.thickness * u"m"
    
    for triangle in DT.each_solid_triangle(tri)
        i, j, k = triangle
        push!(shells, ShellTri3(
            (node_map[i], node_map[j], node_map[k]), 
            thickness_q, 
            section.E * u"Pa", 
            section.ν; 
            ρ=section.ρ,
            id=id
        ))
    end

    # A collinear boundary emits a triangle with no area. Drop it here, and
    # weld corners that are the same point, before any model is assembled.
    collapse_degenerate_shells!(shells)
    return shells
end

"""
Generate mesh points including support line points and interior nodes.
Returns (all_points, support_point_indices, interior_node_point_indices).

The third return value maps input interior node index → output point index,
allowing the caller to reuse actual Node objects instead of creating new ones.
"""
function _generate_mesh_points_with_supports(
    boundary_pts::Vector{Tuple{Float64, Float64}}, 
    n::Int,
    support_lines::Vector{_SupportLine},
    interior_node_pts::Vector{Tuple{Float64, Float64}} = Tuple{Float64, Float64}[],
    constrained_lines::Vector{Vector{Tuple{Float64, Float64}}} = Vector{Tuple{Float64, Float64}}[]
)
    nc = length(boundary_pts)
    
    # Tolerance for deduplication (in meters)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)
    
    # Auto-generated points (grid, edge subdivisions, support samples) that land
    # on a constrained-line segment interior must be dropped, or DT would split
    # the constrained edge there. Boundary corners and interior nodes are exempt.
    _excluded(pt) = !isempty(constrained_lines) &&
        _on_any_constrained_segment_interior(pt, constrained_lines, tol)
    
    # Use a dict to dedupe: key is rounded coords, value is (exact coords, is_support, interior_node_idx)
    # interior_node_idx=0 means not an interior node
    points_dict = Dict{Tuple{Int64, Int64}, Tuple{Tuple{Float64, Float64}, Bool, Int}}()
    
    # Add boundary points first (these take priority)
    for pt in boundary_pts
        key = (round_coord(pt[1]), round_coord(pt[2]))
        points_dict[key] = (pt, false, 0)
    end
    
    # Add interior node points (these are actual nodes we want to embed)
    for (idx, pt) in enumerate(interior_node_pts)
        key = (round_coord(pt[1]), round_coord(pt[2]))
        if !haskey(points_dict, key)
            # Only add if inside polygon or on boundary
            if _point_inside_polygon(pt, boundary_pts) || _is_on_boundary_edge(pt, boundary_pts, tol)
                points_dict[key] = (pt, false, idx)  # Track which interior node
            end
        end
    end
    
    # Add support line points
    for line in support_lines
        line_pts = _points_along_line(line.start_pt, line.end_pt, n)
        for pt in line_pts
            key = (round_coord(pt[1]), round_coord(pt[2]))
            if !haskey(points_dict, key) && !_excluded(pt)
                # Only add if inside polygon
                if _point_inside_polygon(pt, boundary_pts) || _is_on_boundary_edge(pt, boundary_pts, tol)
                    points_dict[key] = (pt, true, 0)  # Mark as support point
                end
            end
        end
    end
    
    if n <= 1
        # Build result with support tracking
        result = collect(boundary_pts)
        support_indices = Set{Int}()
        interior_node_map = Dict{Int, Int}()  # input idx → point idx
        
        for (key, (pt, is_support, interior_idx)) in points_dict
            is_boundary = any(abs(pt[1] - bp[1]) < tol && abs(pt[2] - bp[2]) < tol for bp in boundary_pts)
            if !is_boundary
                push!(result, pt)
                point_idx = length(result)
                if is_support
                    push!(support_indices, point_idx)
                end
                if interior_idx > 0
                    interior_node_map[interior_idx] = point_idx
                end
            end
        end
        
        return result, support_indices, interior_node_map
    end
    
    # Compute bounding box
    xmin = minimum(p[1] for p in boundary_pts)
    xmax = maximum(p[1] for p in boundary_pts)
    ymin = minimum(p[2] for p in boundary_pts)
    ymax = maximum(p[2] for p in boundary_pts)
    
    dx = (xmax - xmin) / n
    dy = (ymax - ymin) / n
    
    # Add edge subdivision points (not supports)
    for i in 1:nc
        p1 = boundary_pts[i]
        p2 = boundary_pts[mod1(i+1, nc)]
        
        for k in 1:n-1
            t = k / n
            pt = (p1[1] * (1-t) + p2[1] * t, p1[2] * (1-t) + p2[2] * t)
            key = (round_coord(pt[1]), round_coord(pt[2]))
            if !haskey(points_dict, key) && !_excluded(pt)
                points_dict[key] = (pt, false, 0)
            end
        end
    end
    
    # Add interior grid points (not supports)
    for j in 1:n-1
        for i in 1:n-1
            x = xmin + i * dx
            y = ymin + j * dy
            pt = (x, y)
            
            if _point_inside_polygon(pt, boundary_pts) && !_excluded(pt)
                key = (round_coord(pt[1]), round_coord(pt[2]))
                if !haskey(points_dict, key)
                    points_dict[key] = (pt, false, 0)
                end
            end
        end
    end
    
    # Build result: boundary points first, then others
    result = collect(boundary_pts)
    support_indices = Set{Int}()
    interior_node_map = Dict{Int, Int}()  # input idx → point idx
    
    for (key, (pt, is_support, interior_idx)) in points_dict
        is_boundary = any(abs(pt[1] - bp[1]) < tol && abs(pt[2] - bp[2]) < tol for bp in boundary_pts)
        if !is_boundary
            push!(result, pt)
            point_idx = length(result)
            if is_support
                push!(support_indices, point_idx)
            end
            if interior_idx > 0
                interior_node_map[interior_idx] = point_idx
            end
        end
    end
    
    return result, support_indices, interior_node_map
end

"""
Generate points along a line at spacing consistent with n divisions.
"""
function _points_along_line(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, n::Int)
    points = Tuple{Float64, Float64}[]
    
    # Include endpoints and n-1 interior points
    for k in 0:n
        t = k / n
        pt = (p1[1] * (1-t) + p2[1] * t, p1[2] * (1-t) + p2[2] * t)
        push!(points, pt)
    end
    
    return points
end

"""Clockwise loop, so a hole boundary has the solid region on its left."""
function _as_hole_loop(raw)
    loop = Tuple{Float64, Float64}[(Float64(p[1]), Float64(p[2])) for p in raw]
    length(loop) > 1 && loop[end] == loop[1] && pop!(loop)
    length(loop) < 3 && return Tuple{Float64, Float64}[]
    area = 0.0
    n = length(loop)
    for i in 1:n
        p = loop[i]
        q = loop[i == n ? 1 : i + 1]
        area += p[1] * q[2] - q[1] * p[2]
    end
    area > 0 && reverse!(loop)
    return loop
end

"""False when `pt` is strictly inside a hole. A point on the hole edge stays."""
function _outside_holes(pt, holes)
    for h in holes
        _point_inside_polygon(pt, h) && !_is_on_boundary_edge(pt, h, 1e-7) && return false
    end
    return true
end

"""Drop mesh points that landed inside a hole. The first `nc` boundary points stay."""
function _drop_hole_interior(points, nc, holes)
    isempty(holes) && return points, Dict(i => i for i in eachindex(points))
    remap = Dict{Int, Int}()
    out = Tuple{Float64, Float64}[]
    for (i, p) in enumerate(points)
        if i <= nc || _outside_holes(p, holes)
            push!(out, p)
            remap[i] = length(out)
        end
    end
    return out, remap
end

"""
Outer boundary plus one clockwise curve per hole.

`DelaunayTriangulation` reads a `Vector{Vector{Int}}` as one curve and a
`Vector{Vector{Vector{Int}}}` as several. Holes are the extra curves.
"""
function _boundary_curves!(points, nc, holes)
    outer = vcat(collect(1:nc), [1])
    isempty(holes) && return Any[outer]
    curves = Vector{Vector{Vector{Int}}}()
    push!(curves, [outer])
    for h in holes
        # Reuse a seed point already on the hole (a wall node, a corner) so
        # the slab and the wall are the same Node. Points that lie on an
        # edge are inserted in order; the rest of the corner is appended.
        push!(curves, [_hole_curve_indices!(points, h)])
    end
    return curves
end

"""Indices of `loop` in `points`, closed, with existing points on each edge inserted."""
function _hole_curve_indices!(points, loop; tol::Float64 = 1e-6)
    function locate(p)
        for (i, q) in enumerate(points)
            hypot(q[1] - p[1], q[2] - p[2]) <= tol && return i
        end
        push!(points, (Float64(p[1]), Float64(p[2])))
        return length(points)
    end
    corners = [locate(p) for p in loop]
    ids = Int[]
    n = length(corners)
    for i in 1:n
        a = corners[i]
        b = corners[mod1(i + 1, n)]
        push!(ids, a)
        pa, pb = points[a], points[b]
        abx, aby = pb[1] - pa[1], pb[2] - pa[2]
        len2 = abx * abx + aby * aby
        len2 == 0 && continue
        mids = Tuple{Float64, Int}[]
        for (k, q) in enumerate(points)
            (k == a || k == b) && continue
            s = ((q[1] - pa[1]) * abx + (q[2] - pa[2]) * aby) / len2
            (0 < s < 1) || continue
            hypot(q[1] - (pa[1] + s * abx), q[2] - (pa[2] + s * aby)) <= tol || continue
            push!(mids, (s, k))
        end
        sort!(mids)
        for (_, k) in mids
            push!(ids, k)
        end
    end
    push!(ids, ids[1])
    return ids
end

_on_hole_edge(pt, holes, tol) = any(h -> _is_on_boundary_edge(pt, h, tol), holes)

function _remap_pairs(pairs, remap)
    out = Set{Tuple{Int, Int}}()
    for (i, j) in pairs
        haskey(remap, i) && haskey(remap, j) && push!(out, (remap[i], remap[j]))
    end
    return out
end

function _point_inside_polygon(pt::Tuple{Float64, Float64}, polygon::Vector{Tuple{Float64, Float64}})
    x, y = pt
    n = length(polygon)
    inside = false
    
    j = n
    for i in 1:n
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    
    return inside
end

# ============================================================================
# Structured mesh → ShellTri3 elements
# ============================================================================

"""
    _create_shells_from_structured(xs, ys, corners, boundary_pts, interior_nodes,
                                    section, id, edge_support_type, interior_support_type)

Build ShellTri3 elements from a structured rectangular grid via `_t3blockx`.
Maps grid nodes to Asap `Node`s with appropriate fixities.
"""
function _create_shells_from_structured(
    xs::Vector{Float64},
    ys::Vector{Float64},
    corners,
    boundary_pts::Vector{Tuple{Float64, Float64}},
    interior_nodes::Vector{Node},
    section::ShellSection,
    id::Symbol,
    edge_support_type::Union{Symbol, Vector{Bool}},
    interior_support_type::Union{Symbol, Vector{Bool}}
)
    # Generate structured mesh
    mesh_pts, mesh_conn = _t3blockx(xs, ys)
    
    # Average z from corners
    nc = length(corners)
    z_val = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
    
    tol = 1e-6
    
    # Build map from interior node positions to actual Node objects (for structural connectivity)
    interior_map = Dict{Tuple{Int64, Int64}, Node}()
    round_coord(x) = round(Int64, x / tol)
    for nd in interior_nodes
        key = (round_coord(ustrip(u"m", nd.position[1])),
               round_coord(ustrip(u"m", nd.position[2])))
        interior_map[key] = nd
    end
    
    # Build map from corner positions to actual corner Node objects
    corner_map = Dict{Tuple{Int64, Int64}, Node}()
    for corner in corners
        key = (round_coord(ustrip(u"m", corner.position[1])),
               round_coord(ustrip(u"m", corner.position[2])))
        corner_map[key] = corner
    end
    
    # Create nodes for each mesh point
    node_array = Vector{Node}(undef, length(mesh_pts))
    
    for (i, pt) in enumerate(mesh_pts)
        key = (round_coord(pt[1]), round_coord(pt[2]))
        
        if haskey(corner_map, key)
            # Use the actual corner node
            node_array[i] = corner_map[key]
        elseif haskey(interior_map, key)
            # Use the actual interior node (structural connectivity)
            node_array[i] = interior_map[key]
        else
            pos = [pt[1], pt[2], z_val] .* u"m"
            
            # Determine fixity
            is_edge = _is_on_boundary_edge(pt, boundary_pts, tol)
            fixity = is_edge ? edge_support_type : :free
            
            node_array[i] = _create_node_with_fixity(pos, fixity, id)
        end
    end
    
    # Create shell elements
    thickness_q = section.thickness * u"m"
    shells = ShellTri3[]
    for (i, j, k) in mesh_conn
        push!(shells, ShellTri3(
            (node_array[i], node_array[j], node_array[k]),
            thickness_q,
            section.E * u"Pa",
            section.ν;
            ρ=section.ρ,
            id=id
        ))
    end
    
    return shells
end

# ============================================================================
# Refinement rings for Delaunay meshing
# ============================================================================

"""
    _add_refinement_rings!(all_points, refine_pts, h_near, r_transition,
                           boundary_pts, patches=ShellPatch[])

Add concentric rings of points around each refinement target to improve
Delaunay mesh density near supports.

Points are discarded if they:
- fall outside the slab boundary polygon,
- duplicate an existing mesh point, or
- fall **inside any ShellPatch** interior (column footprint).  Placing ring
  points inside a patch creates degenerate slivers between the ring nodes
  and the patch's constrained edges, tanking the stiffness-matrix condition
  number.  Patch geometry must be added to `all_points` *before* calling
  this function.
"""
function _add_refinement_rings!(
    all_points::Vector{Tuple{Float64, Float64}},
    refine_pts::Vector{Tuple{Float64, Float64}},
    h_near::Float64,
    r_transition::Float64,
    boundary_pts::Vector{Tuple{Float64, Float64}},
    patches::Vector{ShellPatch} = ShellPatch[]
)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)
    
    # Existing point set for deduplication
    existing = Set{Tuple{Int64, Int64}}()
    for pt in all_points
        push!(existing, (round_coord(pt[1]), round_coord(pt[2])))
    end

    # Precompute patch vertex arrays for point-in-polygon tests
    patch_polys = [p.vertices for p in patches]
    
    for center in refine_pts
        # Add rings at increasing radii
        radii = Float64[]
        r = h_near
        while r < r_transition
            push!(radii, r)
            r *= 1.6  # growth factor
        end
        
        for r in radii
            # Number of points on this ring (circumference / desired spacing)
            n_ring = max(6, round(Int, 2π * r / h_near))
            for k in 0:n_ring-1
                θ = 2π * k / n_ring
                pt = (center[1] + r * cos(θ), center[2] + r * sin(θ))
                key = (round_coord(pt[1]), round_coord(pt[2]))
                
                # Skip duplicates
                key in existing && continue
                # Must be inside the slab boundary
                _point_inside_polygon(pt, boundary_pts) || continue
                # Must NOT be inside any patch interior
                any(poly -> _point_inside_polygon(pt, poly), patch_polys) && continue
                
                push!(all_points, pt)
                push!(existing, key)
            end
        end
    end
end

# ============================================================================
