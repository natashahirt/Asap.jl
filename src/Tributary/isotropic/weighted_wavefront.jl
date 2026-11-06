# =============================================================================
# Weighted Straight Skeleton: Affine Motion Model
# =============================================================================

"""
Given CCW polygon pts[1..n] and per-edge weights w[1..n],
return inward unit normals n[i] and constants c[i] for each edge i.
Edge i is pts[i] -> pts[i+1].
"""
function _edge_lines_ccw(pts::Vector{NTuple{2,Float64}}, w::Vector{Float64})
    n = length(pts)
    normals = Vector{NTuple{2,Float64}}(undef, n)
    cs = Vector{Float64}(undef, n)
    
    for i in 1:n
        p = pts[i]
        q = pts[mod1(i+1, n)]
        ex, ey = (q[1]-p[1], q[2]-p[2])
        el = hypot(ex, ey)
        el < 1e-12 && error("Degenerate edge $i")
        
        ex /= el; ey /= el
        # Inward normal for CCW polygon (90° CCW rotation)
        nx, ny = (-ey, ex)
        normals[i] = (nx, ny)
        cs[i] = nx*p[1] + ny*p[2]
    end
    return normals, cs
end

"""
Given current CCW polygon pts and inward unit normals normals[i] for edges,
compute vertex velocities v[i] satisfying:
n_{i-1}·v[i] = w_{i-1}
n_i·v[i]     = w_i

Note: x0 is the current vertex positions (t=0 state), so we only compute v.
"""
function _vertex_affine_velocity(normals::Vector{NTuple{2,Float64}}, w::Vector{Float64}; det_tol=1e-12)
    n = length(normals)
    v = Vector{NTuple{2,Float64}}(undef, n)
    
    for i in 1:n
        im1 = mod1(i-1, n)
        n1x, n1y = normals[im1]
        n2x, n2y = normals[i]
        
        detA = n1x*n2y - n1y*n2x
        if abs(detA) < det_tol
            # Don't error here; mark as degenerate and let cleanup remove it
            v[i] = (NaN, NaN)
            continue
        end
        
        invA11 =  n2y/detA
        invA12 = -n1y/detA
        invA21 = -n2x/detA
        invA22 =  n1x/detA
        
        w1 = w[im1]
        w2 = w[i]
        vx = invA11*w1 + invA12*w2
        vy = invA21*w1 + invA22*w2
        v[i] = (vx, vy)
    end
    return v
end

