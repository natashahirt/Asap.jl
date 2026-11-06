# =============================================================================
# Wavefront Cleanup (removes degenerate edges/vertices)
# =============================================================================

"""Cross product helper: cross of (b-a) x (c-b)."""
function _cross(a, b, c)
    return (b[1]-a[1])*(c[2]-b[2]) - (b[2]-a[2])*(c[1]-b[1])
end

"""
Clean wavefront polygon to avoid degenerate edges/vertices.

Removes:
1. Consecutive duplicate points / tiny edges
2. Near-collinear vertices (180° turns)
3. Tiny edges introduced by collinear removal

Returns (pts, edge_map, weight_map) with consistent lengths.
Assumes pts are CCW and represent a simple cycle.
"""
function cleanup_wavefront(pts::Vector{NTuple{2,Float64}},
                           edge_map::Vector{Int},
                           w::Vector{Float64};
                           eps_len=1e-10,
                           eps_col=1e-12)
    n = length(pts)
    n == length(edge_map) == length(w) || error("cleanup: length mismatch")
    n < 3 && return pts, edge_map, w
    
    # 1) Remove consecutive duplicates / tiny edges
    keep = trues(n)
    for i in 1:n
        j = mod1(i+1, n)
        if _dist(pts[i], pts[j]) < eps_len
            # Remove vertex j (equivalently remove edge i)
            keep[j] = false
        end
    end
    
    pts2 = NTuple{2,Float64}[]
    em2 = Int[]
    w2 = Float64[]
    for i in 1:n
        if keep[i]
            push!(pts2, pts[i])
            push!(em2, edge_map[i])
            push!(w2, w[i])
        end
    end
    
    # After dropping vertices, ensure cycle consistency
    n = length(pts2)
    n < 3 && return pts2, em2, w2
    
    # 2) Remove near-collinear vertices (180° turns)
    # Vertex i is between edges (i-1) and i.
    keep = trues(n)
    for i in 1:n
        a = pts2[mod1(i-1, n)]
        b = pts2[i]
        c = pts2[mod1(i+1, n)]
        if abs(_cross(a, b, c)) < eps_col
            keep[i] = false
        end
    end
    
    pts3 = NTuple{2,Float64}[]
    em3 = Int[]
    w3 = Float64[]
    for i in 1:n
        if keep[i]
            push!(pts3, pts2[i])
            push!(em3, em2[i])
            push!(w3, w2[i])
        end
    end
    
    n = length(pts3)
    n < 3 && return pts3, em3, w3
    
    # 3) One more pass removing tiny edges introduced by collinear removal
    keep = trues(n)
    for i in 1:n
        j = mod1(i+1, n)
        if _dist(pts3[i], pts3[j]) < eps_len
            keep[j] = false
        end
    end
    
    pts4 = NTuple{2,Float64}[]
    em4 = Int[]
    w4 = Float64[]
    for i in 1:n
        if keep[i]
            push!(pts4, pts3[i])
            push!(em4, em3[i])
            push!(w4, w3[i])
        end
    end
    
    return pts4, em4, w4
end

