# ============================================================================
# Polygon Validation Helpers
# ============================================================================

"""Shoelace formula for signed polygon area (positive if CCW)."""
function _shoelace_area(pts::Vector{<:Tuple{<:Real, <:Real}})
    n = length(pts)
    n < 3 && return 0.0
    area = sum(pts[i][1] * pts[mod1(i+1,n)][2] - pts[mod1(i+1,n)][1] * pts[i][2] for i in 1:n) / 2
    return area
end

"""Check if two line segments intersect (not at endpoints)."""
function _segments_intersect(
    p1::Tuple{<:Real, <:Real}, p2::Tuple{<:Real, <:Real},
    p3::Tuple{<:Real, <:Real}, p4::Tuple{<:Real, <:Real}
)
    # Direction vectors
    d1 = (p2[1] - p1[1], p2[2] - p1[2])
    d2 = (p4[1] - p3[1], p4[2] - p3[2])
    
    # Cross product (2D)
    cross = d1[1] * d2[2] - d1[2] * d2[1]
    abs(cross) < 1e-12 && return false  # Parallel
    
    # Solve for parameters t and s
    t = ((p3[1] - p1[1]) * d2[2] - (p3[2] - p1[2]) * d2[1]) / cross
    s = ((p3[1] - p1[1]) * d1[2] - (p3[2] - p1[2]) * d1[1]) / cross
    
    # Intersection in interior of both segments (exclude endpoints)
    ε = 1e-9
    return ε < t < 1-ε && ε < s < 1-ε
end

"""Check if polygon is self-intersecting (edges cross each other)."""
function _is_self_intersecting(pts::Vector{<:Tuple{<:Real, <:Real}})
    n = length(pts)
    n < 4 && return false  # Triangle can't self-intersect
    
    for i in 1:n
        p1, p2 = pts[i], pts[mod1(i+1, n)]
        # Check against non-adjacent edges
        for j in (i+2):n
            j == n && i == 1 && continue  # Skip edge that shares vertex with edge i
            p3, p4 = pts[j], pts[mod1(j+1, n)]
            _segments_intersect(p1, p2, p3, p4) && return true
        end
    end
    return false
end

