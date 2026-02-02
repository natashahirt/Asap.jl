#=
Shell-to-Beam Load Distribution (Internal Implementation)
=========================================================

Internal functions for converting shell panel geometry into TributaryLoads.
These are implementation details used by AreaLoad when distribute_to=elements.

Public API: Use `AreaLoad(shells, pressure; distribute_to=beams)` instead.

Attribution: Part of the ASAP structural analysis framework.
=#

using Meshes: Point, coords
using Unitful: ustrip, @u_str, Quantity

# =============================================================================
# Edge Matching (Internal)
# =============================================================================

"""Internal: Match beam to polygon edge index."""
function _match_beam_to_polygon_edge(
    beam::FrameElement,
    vertices::Vector{<:Point};
    tol::Float64 = 1e-6
)::Union{Int, Nothing}
    n = length(vertices)
    n < 2 && return nothing
    
    # Get beam endpoint positions in meters (Node stores position as a vector)
    # Element uses nodeStart and nodeEnd fields
    beam_start = (
        Float64(ustrip(u"m", beam.nodeStart.position[1])),
        Float64(ustrip(u"m", beam.nodeStart.position[2]))
    )
    beam_end = (
        Float64(ustrip(u"m", beam.nodeEnd.position[1])),
        Float64(ustrip(u"m", beam.nodeEnd.position[2]))
    )
    
    # Check each polygon edge
    for i in 1:n
        j = mod1(i + 1, n)
        
        # Get edge endpoint positions
        c_i = coords(vertices[i])
        c_j = coords(vertices[j])
        v1 = (Float64(ustrip(u"m", c_i.x)), Float64(ustrip(u"m", c_i.y)))
        v2 = (Float64(ustrip(u"m", c_j.x)), Float64(ustrip(u"m", c_j.y)))
        
        # Check both orientations
        dist_forward = hypot(beam_start[1] - v1[1], beam_start[2] - v1[2]) +
                       hypot(beam_end[1] - v2[1], beam_end[2] - v2[2])
        dist_reverse = hypot(beam_start[1] - v2[1], beam_start[2] - v2[2]) +
                       hypot(beam_end[1] - v1[1], beam_end[2] - v1[2])
        
        if dist_forward < 2 * tol || dist_reverse < 2 * tol
            return i
        end
    end
    
    return nothing
end

"""Internal: Match beam (as tuple coords) to polygon edge index."""
function _match_beam_to_polygon_edge(
    beam_start::NTuple{2, Float64},
    beam_end::NTuple{2, Float64},
    vertices::Vector{NTuple{2, Float64}};
    tol::Float64 = 1e-6
)::Union{Int, Nothing}
    n = length(vertices)
    n < 2 && return nothing
    
    for i in 1:n
        j = mod1(i + 1, n)
        v1, v2 = vertices[i], vertices[j]
        
        dist_forward = hypot(beam_start[1] - v1[1], beam_start[2] - v1[2]) +
                       hypot(beam_end[1] - v2[1], beam_end[2] - v2[2])
        dist_reverse = hypot(beam_start[1] - v2[1], beam_start[2] - v2[2]) +
                       hypot(beam_end[1] - v1[1], beam_end[2] - v1[2])
        
        if dist_forward < 2 * tol || dist_reverse < 2 * tol
            return i
        end
    end
    
    return nothing
end

# =============================================================================
# Shell Panel Boundary Extraction (Internal)
# =============================================================================

"""Internal: Extract outer boundary polygon from shell elements."""
function _get_shell_panel_boundary(
    shells::Vector{<:ShellElement};
    panel_id::Union{Symbol, Nothing} = nothing
)::Vector{NTuple{2, Float64}}
    # Filter by panel_id if specified
    filtered = isnothing(panel_id) ? shells : filter(s -> s.id == panel_id, shells)
    isempty(filtered) && return NTuple{2, Float64}[]
    
    # Collect all edges (as sorted node position pairs to identify shared edges)
    # Edge format: ((x1,y1), (x2,y2)) with (x1,y1) < (x2,y2) lexicographically
    edge_count = Dict{Tuple{NTuple{2,Float64}, NTuple{2,Float64}}, Int}()
    edge_to_original = Dict{Tuple{NTuple{2,Float64}, NTuple{2,Float64}}, NTuple{2, NTuple{2,Float64}}}()
    
    for shell in filtered
        nodes = shell.nodes
        n_nodes = length(nodes)
        
        for i in 1:n_nodes
            j = mod1(i + 1, n_nodes)
            n1, n2 = nodes[i], nodes[j]
            
            # Node stores position as a vector: position[1]=x, position[2]=y, position[3]=z
            p1 = (Float64(ustrip(u"m", n1.position[1])), Float64(ustrip(u"m", n1.position[2])))
            p2 = (Float64(ustrip(u"m", n2.position[1])), Float64(ustrip(u"m", n2.position[2])))
            
            # Create canonical edge key (sorted)
            key = p1 < p2 ? (p1, p2) : (p2, p1)
            original = (p1, p2)  # Keep original orientation
            
            edge_count[key] = get(edge_count, key, 0) + 1
            if !haskey(edge_to_original, key)
                edge_to_original[key] = original
            end
        end
    end
    
    # Keep only boundary edges (count == 1)
    boundary_edges = [(edge_to_original[k], k) for (k, count) in edge_count if count == 1]
    isempty(boundary_edges) && return NTuple{2, Float64}[]
    
    # Chain edges into a polygon
    # Build adjacency: point → list of edges starting/ending there
    point_to_edges = Dict{NTuple{2,Float64}, Vector{Int}}()
    edges = [e[1] for e in boundary_edges]  # Original orientation
    
    for (idx, (p1, p2)) in enumerate(edges)
        push!(get!(point_to_edges, p1, Int[]), idx)
        push!(get!(point_to_edges, p2, Int[]), idx)
    end
    
    # Walk the boundary starting from first edge
    polygon = NTuple{2, Float64}[]
    used = Set{Int}()
    
    current_edge_idx = 1
    current_point = edges[1][1]
    push!(polygon, current_point)
    push!(used, 1)
    current_point = edges[1][2]
    push!(polygon, current_point)
    
    while length(used) < length(edges)
        # Find next edge that shares current_point
        found = false
        for edge_idx in get(point_to_edges, current_point, Int[])
            if edge_idx ∉ used
                push!(used, edge_idx)
                p1, p2 = edges[edge_idx]
                # Move to the other end of this edge
                next_point = (p1 == current_point) ? p2 : p1
                push!(polygon, next_point)
                current_point = next_point
                found = true
                break
            end
        end
        !found && break
    end
    
    # Remove duplicate last point if it matches first
    if length(polygon) > 1 && polygon[end] == polygon[1]
        pop!(polygon)
    end
    
    # Ensure CCW orientation
    if _polygon_area(polygon) < 0
        reverse!(polygon)
    end
    
    return polygon
end

# Public wrappers for backward compatibility / testing
"""Match beam to polygon edge. See AreaLoad for preferred API."""
match_beam_to_polygon_edge(beam::FrameElement, vertices::Vector{<:Point}; tol=1e-6) = 
    _match_beam_to_polygon_edge(beam, vertices; tol=tol)
match_beam_to_polygon_edge(start::NTuple{2,Float64}, stop::NTuple{2,Float64}, verts::Vector{NTuple{2,Float64}}; tol=1e-6) = 
    _match_beam_to_polygon_edge(start, stop, verts; tol=tol)

"""Extract shell panel boundary. See AreaLoad for preferred API."""
get_shell_panel_boundary(shells::Vector{<:ShellElement}; panel_id=nothing) = 
    _get_shell_panel_boundary(shells; panel_id=panel_id)

# =============================================================================
# Helper: Extract Depth Profile from Tributary Polygon
# =============================================================================

"""Internal: Extract piecewise-linear depth profile from tributary polygon vertices."""
function _extract_depth_profile(tp)
    # Empty polygon = zero tributary width
    if isempty(tp.s) || isempty(tp.d)
        return [0.0, 1.0], [0.0, 0.0]
    end
    
    # Sort ALL polygon vertices by s to get the depth profile
    # This works because the polygon boundary traces from beam edge (d=0) 
    # up to interior edge (d>0) and back down
    perm = sortperm(tp.s)
    s_sorted = tp.s[perm]
    d_sorted = tp.d[perm]
    
    # Remove duplicate (s, d) pairs 
    positions = Float64[]
    depths = Float64[]
    for i in eachindex(s_sorted)
        # Take the MAXIMUM d at each s (upper envelope)
        if isempty(positions) || s_sorted[i] > positions[end] + 1e-9
            push!(positions, s_sorted[i])
            push!(depths, max(0.0, d_sorted[i]))
        else
            # Same s position - take max d
            depths[end] = max(depths[end], d_sorted[i])
        end
    end
    
    # Ensure we have endpoints at s=0 and s=1
    if isempty(positions) || positions[1] > 0.01
        pushfirst!(positions, 0.0)
        pushfirst!(depths, 0.0)
    end
    if positions[end] < 0.99
        push!(positions, 1.0)
        push!(depths, 0.0)
    end
    
    return positions, depths
end

# =============================================================================
# Shell to Tributary Load Conversion (Internal)
# =============================================================================

"""Internal: Main implementation for shell-to-tributary conversion."""
function _shell_to_tributary_loads(
    shells::Vector{<:ShellElement},
    beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
    tol::Float64 = 1e-6,
    panel_id::Union{Symbol, Nothing} = nothing,
    interior_beams::Vector{<:FrameElement} = FrameElement[]
)::Vector{TributaryLoad}
    isempty(shells) && return TributaryLoad[]
    isempty(beams) && return TributaryLoad[]
    
    # Extract panel boundary
    boundary = _get_shell_panel_boundary(shells; panel_id=panel_id)
    isempty(boundary) && error("Could not extract panel boundary from shells")
    
    # If interior beams specified, split polygon and compute sub-panel tributaries
    if !isempty(interior_beams)
        return _shell_to_tributary_loads_with_interior(
            boundary, beams, interior_beams, pressure;
            axis=axis, direction=direction, tol=tol
        )
    end
    
    # Convert boundary to Meshes.Point for tributary computation
    boundary_points = [Point(p[1]*u"m", p[2]*u"m") for p in boundary]
    
    # Compute tributary polygons
    trib_polys = isnothing(axis) ? 
        get_tributary_polygons(boundary_points) :
        get_tributary_polygons(boundary_points; axis=axis)
    
    # Map each beam to its edge
    beam_to_edge = Dict{FrameElement, Int}()
    unmatched_beams = FrameElement[]
    
    for beam in beams
        edge_idx = _match_beam_to_polygon_edge(beam, boundary_points; tol=tol)
        if isnothing(edge_idx)
            push!(unmatched_beams, beam)
        else
            beam_to_edge[beam] = edge_idx
        end
    end
    
    # Error if any beam doesn't match
    if !isempty(unmatched_beams)
        beam_ids = [b.id for b in unmatched_beams]
        error("Beams do not match any panel edge: $beam_ids. " *
              "Ensure all beams lie on the shell panel boundary, " *
              "or specify interior_beams for beams inside the panel.")
    end
    
    # Create TributaryLoad for each beam
    loads = TributaryLoad[]
    
    for (beam, edge_idx) in beam_to_edge
        tp = trib_polys[edge_idx]
        
        # Extract depth profile from polygon vertices
        positions, depths = _extract_depth_profile(tp)
        
        # Clamp negative depths to zero (can occur with numerical precision issues)
        depths = max.(depths, 0.0)
        widths = [w * u"m" for w in depths]
        
        trib_load = TributaryLoad(
            beam,
            positions,
            widths,
            pressure,
            direction
        )
        push!(loads, trib_load)
    end
    
    return loads
end

"""
Internal: Handle tributary computation with interior beams.
Splits the polygon along interior beam lines and computes tributary for each sub-panel.
Interior beams receive load from both sides.
"""
function _shell_to_tributary_loads_with_interior(
    boundary::Vector{NTuple{2, Float64}},
    beams::Vector{<:FrameElement},
    interior_beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
    tol::Float64 = 1e-6
)::Vector{TributaryLoad}
    
    # Split the boundary polygon along interior beam lines
    sub_polygons = _split_polygon_by_beams(boundary, interior_beams; tol=tol)
    
    # For each sub-polygon, compute tributary to edge and interior beams
    # Interior beams will accumulate load from adjacent sub-polygons
    beam_loads = Dict{FrameElement, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}()
    
    for sub_poly in sub_polygons
        boundary_points = [Point(p[1]*u"m", p[2]*u"m") for p in sub_poly]
        
        # Compute tributary for this sub-panel
        trib_polys = isnothing(axis) ? 
            get_tributary_polygons(boundary_points) :
            get_tributary_polygons(boundary_points; axis=axis)
        
        # Match beams to sub-panel edges
        for beam in beams
            edge_idx = _match_beam_to_polygon_edge(beam, boundary_points; tol=tol)
            if !isnothing(edge_idx)
                tp = trib_polys[edge_idx]
                positions, depths = _extract_depth_profile(tp)
                depths = max.(depths, 0.0)
                
                if !haskey(beam_loads, beam)
                    beam_loads[beam] = Tuple{Vector{Float64}, Vector{Float64}}[]
                end
                push!(beam_loads[beam], (positions, depths))
            end
        end
    end
    
    # Create TributaryLoads, combining contributions for interior beams
    loads = TributaryLoad[]
    
    for (beam, contributions) in beam_loads
        if length(contributions) == 1
            # Edge beam - single contribution
            positions, depths = contributions[1]
            widths = [w * u"m" for w in depths]
            push!(loads, TributaryLoad(beam, positions, widths, pressure, direction))
        else
            # Interior beam - combine contributions from both sides
            combined = _combine_tributary_depths(contributions)
            widths = [w * u"m" for w in combined[2]]
            push!(loads, TributaryLoad(beam, combined[1], widths, pressure, direction))
        end
    end
    
    return loads
end

"""
Internal: Split a polygon along interior beam lines.
Returns a vector of sub-polygons.
"""
function _split_polygon_by_beams(
    polygon::Vector{NTuple{2, Float64}},
    beams::Vector{<:FrameElement};
    tol::Float64 = 1e-6
)::Vector{Vector{NTuple{2, Float64}}}
    
    # Start with the original polygon
    current_polygons = [polygon]
    
    for beam in beams
        # Get beam endpoints
        p1 = (Float64(ustrip(u"m", beam.nodeStart.position[1])),
              Float64(ustrip(u"m", beam.nodeStart.position[2])))
        p2 = (Float64(ustrip(u"m", beam.nodeEnd.position[1])),
              Float64(ustrip(u"m", beam.nodeEnd.position[2])))
        
        # Split each current polygon by this beam line
        new_polygons = Vector{NTuple{2, Float64}}[]
        
        for poly in current_polygons
            # Try to split this polygon
            split_result = _split_polygon_by_line(poly, p1, p2; tol=tol)
            append!(new_polygons, split_result)
        end
        
        current_polygons = new_polygons
    end
    
    return current_polygons
end

"""
Internal: Split a polygon by a line segment.
If the line doesn't intersect the polygon interior, returns [polygon].
If it does, returns the two resulting sub-polygons.
"""
function _split_polygon_by_line(
    polygon::Vector{NTuple{2, Float64}},
    line_start::NTuple{2, Float64},
    line_end::NTuple{2, Float64};
    tol::Float64 = 1e-6
)::Vector{Vector{NTuple{2, Float64}}}
    
    n = length(polygon)
    n < 3 && return [polygon]
    
    # Find intersections of the line with polygon edges
    intersections = Tuple{Int, Float64, NTuple{2, Float64}}[]  # (edge_idx, t_along_edge, point)
    
    for i in 1:n
        j = mod1(i + 1, n)
        e1, e2 = polygon[i], polygon[j]
        
        result = _line_segment_intersection(line_start, line_end, e1, e2; tol=tol)
        if !isnothing(result)
            t, pt = result
            push!(intersections, (i, t, pt))
        end
    end
    
    # Need exactly 2 intersections to split
    length(intersections) != 2 && return [polygon]
    
    # Sort by edge index
    sort!(intersections, by=x->x[1])
    
    idx1, t1, pt1 = intersections[1]
    idx2, t2, pt2 = intersections[2]
    
    # Build two sub-polygons
    # Polygon A: from pt1 along edges to pt2, then back via the split line
    # Polygon B: from pt2 along edges to pt1, then back via the split line
    
    poly_a = NTuple{2, Float64}[]
    poly_b = NTuple{2, Float64}[]
    
    # Add first intersection point to both
    push!(poly_a, pt1)
    push!(poly_b, pt1)
    
    # Walk from edge idx1+1 to idx2
    i = mod1(idx1 + 1, n)
    while i != mod1(idx2 + 1, n)
        push!(poly_a, polygon[i])
        i = mod1(i + 1, n)
    end
    
    # Add second intersection point to poly_a, then close with pt1
    push!(poly_a, pt2)
    
    # Walk from edge idx2+1 back to idx1
    i = mod1(idx2 + 1, n)
    while i != mod1(idx1 + 1, n)
        push!(poly_b, polygon[i])
        i = mod1(i + 1, n)
    end
    
    # Add pt2 to poly_b to close via split line
    push!(poly_b, pt2)
    
    # Return both if they have enough vertices
    result = Vector{NTuple{2, Float64}}[]
    length(poly_a) >= 3 && push!(result, poly_a)
    length(poly_b) >= 3 && push!(result, poly_b)
    
    isempty(result) && return [polygon]
    return result
end

"""
Internal: Find intersection of two line segments.
Returns (t, point) where t is parameter along first segment, or nothing if no intersection.
"""
function _line_segment_intersection(
    a1::NTuple{2, Float64}, a2::NTuple{2, Float64},
    b1::NTuple{2, Float64}, b2::NTuple{2, Float64};
    tol::Float64 = 1e-6
)::Union{Nothing, Tuple{Float64, NTuple{2, Float64}}}
    
    dx_a = a2[1] - a1[1]
    dy_a = a2[2] - a1[2]
    dx_b = b2[1] - b1[1]
    dy_b = b2[2] - b1[2]
    
    denom = dx_a * dy_b - dy_a * dx_b
    abs(denom) < tol && return nothing  # Parallel
    
    dx_ab = b1[1] - a1[1]
    dy_ab = b1[2] - a1[2]
    
    t = (dx_ab * dy_b - dy_ab * dx_b) / denom
    s = (dx_ab * dy_a - dy_ab * dx_a) / denom
    
    # Check if intersection is within both segments (with tolerance at endpoints)
    if t >= -tol && t <= 1 + tol && s >= -tol && s <= 1 + tol
        t_clamped = clamp(t, 0.0, 1.0)
        pt = (a1[1] + t_clamped * dx_a, a1[2] + t_clamped * dy_a)
        return (s, pt)  # Return s (position along polygon edge) and intersection point
    end
    
    return nothing
end

# =============================================================================
# Public API (Deprecated - use AreaLoad instead)
# =============================================================================

"""
    shell_to_tributary_loads(shells, beams, pressure; axis=nothing, ...)

DEPRECATED: Use `AreaLoad(shells, pressure; distribute_to=beams)` instead.

Convert shell panel geometry into TributaryLoads for edge beams.
"""
function shell_to_tributary_loads(
    shells::Vector{<:ShellElement},
    beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    kwargs...
)
    _shell_to_tributary_loads(shells, beams, pressure; axis=axis, kwargs...)
end

function shell_to_tributary_loads(
    shells::Vector{<:ShellElement},
    beams::Vector{<:FrameElement},
    pressure::Quantity,
    axis::NTuple{2, <:Real};
    kwargs...
)
    axis_vec = [Float64(axis[1]), Float64(axis[2])]
    _shell_to_tributary_loads(shells, beams, pressure; axis=axis_vec, kwargs...)
end

function shell_to_tributary_loads(
    shells::Vector{<:ShellElement},
    beams::Vector{<:FrameElement},
    pressure::Quantity,
    axis_symbol::Symbol;
    kwargs...
)
    axis = if axis_symbol == :x
        [1.0, 0.0]
    elseif axis_symbol == :y
        [0.0, 1.0]
    elseif axis_symbol == :isotropic
        nothing
    else
        error("Unknown axis symbol: $axis_symbol. Use tuple (x, y) or :isotropic")
    end
    _shell_to_tributary_loads(shells, beams, pressure; axis=axis, kwargs...)
end

"""
    shell_panels_to_tributary_loads(panels, beams, pressure; ...)

DEPRECATED: Use `AreaLoad(panels, pressure; distribute_to=beams)` instead.

Generate TributaryLoads from multiple shell panels.
"""
function _shell_panels_to_tributary_loads(
    panels::Vector{Vector{S}},
    beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
    tol::Float64 = 1e-6
)::Vector{TributaryLoad} where S <: ShellElement
    isempty(panels) && return TributaryLoad[]
    isempty(beams) && return TributaryLoad[]
    
    # Track which beams are matched by at least one panel
    beam_matched = Dict{FrameElement, Bool}(b => false for b in beams)
    
    # Accumulate tributary data per beam: beam → list of (positions, widths) tuples
    beam_tributaries = Dict{FrameElement, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}()
    
    for panel_shells in panels
        isempty(panel_shells) && continue
        
        # Extract panel boundary
        boundary = _get_shell_panel_boundary(panel_shells)
        isempty(boundary) && continue
        
        boundary_points = [Point(p[1]*u"m", p[2]*u"m") for p in boundary]
        
        # Compute tributary polygons for this panel
        trib_polys = isnothing(axis) ? 
            get_tributary_polygons(boundary_points) :
            get_tributary_polygons(boundary_points; axis=axis)
        
        # Match beams to this panel's edges
        for beam in beams
            edge_idx = _match_beam_to_polygon_edge(beam, boundary_points; tol=tol)
            if !isnothing(edge_idx)
                beam_matched[beam] = true
                tp = trib_polys[edge_idx]
                
                # Extract depth profile from polygon vertices
                positions, depths = _extract_depth_profile(tp)
                
                if !haskey(beam_tributaries, beam)
                    beam_tributaries[beam] = Tuple{Vector{Float64}, Vector{Float64}}[]
                end
                push!(beam_tributaries[beam], (positions, depths))
            end
        end
    end
    
    # Check for unmatched beams
    unmatched = [b for (b, matched) in beam_matched if !matched]
    if !isempty(unmatched)
        beam_ids = [b.id for b in unmatched]
        error("Beams do not match any panel edge: $beam_ids. " *
              "All beams must lie on at least one panel boundary.")
    end
    
    # Create combined TributaryLoads
    loads = TributaryLoad[]
    
    for beam in beams
        tributaries = get(beam_tributaries, beam, Tuple{Vector{Float64}, Vector{Float64}}[])
        
        if isempty(tributaries)
            # Beam matched but got zero tributary (shouldn't happen after the check above)
            continue
        elseif length(tributaries) == 1
            # Single panel contribution
            positions, depths = tributaries[1]
        else
            # Multiple panels - combine by summing depths at each position
            # Use union of all positions, interpolate depths
            combined = _combine_tributary_depths(tributaries)
            positions, depths = combined
        end
        
        # Clamp negative depths to zero (can occur with numerical precision issues)
        depths = max.(depths, 0.0)
        widths = [d * u"m" for d in depths]
        
        # Handle zero-load case (empty or all-zero depths)
        max_depth = isempty(depths) ? 0.0 : maximum(depths)
        if max_depth < 1e-9
            positions = [0.0, 1.0]
            widths = [0.0u"m", 0.0u"m"]
        end
        
        trib_load = TributaryLoad(
            beam,
            positions,
            widths,
            pressure,
            direction
        )
        push!(loads, trib_load)
    end
    
    return loads
end

# Public wrapper for backward compatibility
function shell_panels_to_tributary_loads(
    panels::Vector{Vector{S}},
    beams::Vector{<:FrameElement},
    pressure::Quantity;
    kwargs...
) where S <: ShellElement
    _shell_panels_to_tributary_loads(panels, beams, pressure; kwargs...)
end

"""Internal: Combine multiple tributary depth profiles by summing depths at each position."""
function _combine_tributary_depths(
    tributaries::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
)::Tuple{Vector{Float64}, Vector{Float64}}
    # Collect all unique positions
    all_positions = Set{Float64}()
    for (positions, _) in tributaries
        union!(all_positions, positions)
    end
    
    sorted_positions = sort(collect(all_positions))
    
    # Sum depths at each position (interpolate if needed)
    combined_depths = zeros(length(sorted_positions))
    
    for (positions, depths) in tributaries
        for (i, s) in enumerate(sorted_positions)
            combined_depths[i] += _interpolate_depth(s, positions, depths)
        end
    end
    
    return (sorted_positions, combined_depths)
end

"""
Linear interpolation of depth at position s.
"""
function _interpolate_depth(s::Float64, positions::Vector{Float64}, depths::Vector{Float64})::Float64
    n = length(positions)
    n == 0 && return 0.0
    n == 1 && return depths[1]
    
    # Clamp to range
    s <= positions[1] && return depths[1]
    s >= positions[end] && return depths[end]
    
    # Find bracketing interval
    for i in 1:(n-1)
        if positions[i] <= s <= positions[i+1]
            t = (s - positions[i]) / (positions[i+1] - positions[i])
            return depths[i] + t * (depths[i+1] - depths[i])
        end
    end
    
    return 0.0
end
