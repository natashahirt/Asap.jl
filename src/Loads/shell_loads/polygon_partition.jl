
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
    
    # The other piece walks the same vertices, but from pt1 the boundary
    # reaches them in reverse. Forward order bows the piece across itself
    # (a trapezoid cut by a joist came back as a bowtie and the one-way
    # partition then counted its area twice).
    mid = NTuple{2, Float64}[]
    i = mod1(idx2 + 1, n)
    while i != mod1(idx1 + 1, n)
        push!(mid, polygon[i])
        i = mod1(i + 1, n)
    end
    append!(poly_b, reverse(mid))
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
