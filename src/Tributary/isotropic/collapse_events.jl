"""
Return t where (x0a + va*t) == (x0b + vb*t).
Least-squares scalar t with residual check.
"""
function _meet_time_affine(x0a, va, x0b, vb; tol_rel=1e-10, tol_abs=1e-12)
    # Check for NaN inputs
    if !isfinite(va[1]) || !isfinite(va[2]) || !isfinite(vb[1]) || !isfinite(vb[2])
        return Inf
    end
    
    ax = va[1] - vb[1]
    ay = va[2] - vb[2]
    bx = x0b[1] - x0a[1]
    by = x0b[2] - x0a[2]
    
    denom = ax*ax + ay*ay
    denom < tol_abs && return Inf
    
    t = (bx*ax + by*ay)/denom
    t <= 1e-12 && return Inf
    
    rx = (x0a[1] + va[1]*t) - (x0b[1] + vb[1]*t)
    ry = (x0a[2] + va[2]*t) - (x0b[2] + vb[2]*t)
    r = hypot(rx, ry)
    
    scale = max(hypot(bx,by), hypot(va[1]*t, va[2]*t), hypot(vb[1]*t, vb[2]*t), 1.0)
    tol = max(tol_abs, tol_rel*scale)
    
    return (r <= tol) ? t : Inf
end

"""
Find all edge collapse events using affine motion model.
Returns (t_min, collapses) where collapses is list of (edge_idx, collapse_point).
"""
function _find_all_collapses_affine(x0::Vector{NTuple{2,Float64}}, v::Vector{NTuple{2,Float64}}, n::Int)
    # First pass: find t_min (skip degenerate vertices with NaN velocities)
    t_min = Inf
    for i in 1:n
        next_i = mod1(i + 1, n)
        # Skip if either vertex has NaN velocity
        if !isfinite(v[i][1]) || !isfinite(v[i][2]) || !isfinite(v[next_i][1]) || !isfinite(v[next_i][2])
            continue
        end
        t = _meet_time_affine(x0[i], v[i], x0[next_i], v[next_i])
        if isfinite(t) && t > 1e-10 && t < t_min
            t_min = t
        end
    end
    
    t_min == Inf && return (Inf, Tuple{Int, NTuple{2,Float64}}[])
    
    # Second pass: collect all edges collapsing at t_min
    tol = max(1e-9, 1e-6 * t_min)
    collapses = Tuple{Int, NTuple{2,Float64}}[]
    
    for i in 1:n
        next_i = mod1(i + 1, n)
        # Skip if either vertex has NaN velocity
        if !isfinite(v[i][1]) || !isfinite(v[i][2]) || !isfinite(v[next_i][1]) || !isfinite(v[next_i][2])
            continue
        end
        t = _meet_time_affine(x0[i], v[i], x0[next_i], v[next_i])
        
        if isfinite(t) && t > 1e-10 && abs(t - t_min) <= tol
            pt = (x0[i][1] + v[i][1]*t, x0[i][2] + v[i][2]*t)
            push!(collapses, (i, pt))
        end
    end
    
    return (t_min, collapses)
end

"""Find all edges collapsing at the minimum time."""
function _find_all_collapses(pts::Vector{NTuple{2,Float64}}, n::Int, bisectors, speeds)
    # First pass: find t_min
    t_min = Inf
    for i in 1:n
        next_i = mod1(i + 1, n)
        p1 = pts[i]
        p2 = pts[next_i]
        d1 = bisectors[i] .* speeds[i]
        d2 = bisectors[next_i] .* speeds[next_i]
        t = _ray_ray_intersect_time(p1, d1, p2, d2)
        if t > 1e-10 && t < t_min
            t_min = t
        end
    end
    
    t_min == Inf && return (Inf, Tuple{Int, NTuple{2,Float64}}[])
    
    # Second pass: collect all edges collapsing at t_min
    tol = max(1e-9, 1e-6 * t_min)
    collapses = Tuple{Int, NTuple{2,Float64}}[]
    
    for i in 1:n
        next_i = mod1(i + 1, n)
        p1 = pts[i]
        p2 = pts[next_i]
        d1 = bisectors[i] .* speeds[i]
        d2 = bisectors[next_i] .* speeds[next_i]
        t = _ray_ray_intersect_time(p1, d1, p2, d2)
        
        if isfinite(t) && t > 1e-10 && abs(t - t_min) <= tol
            pt = (p1[1] + d1[1] * t, p1[2] + d1[2] * t)
            push!(collapses, (i, pt))
        end
    end
    
    return (t_min, collapses)
end

"""Find next single collapse (for final convergence)."""
function _find_next_collapse_single(pts::Vector{NTuple{2,Float64}}, n::Int, bisectors, speeds)
    t_min = Inf
    collapse_idx = 0
    collapse_pt = (0.0, 0.0)
    
    for i in 1:n
        next_i = mod1(i + 1, n)
        p1 = pts[i]
        p2 = pts[next_i]
        d1 = bisectors[i] .* speeds[i]
        d2 = bisectors[next_i] .* speeds[next_i]
        t = _ray_ray_intersect_time(p1, d1, p2, d2)
        
        if t > 1e-10 && t < t_min
            t_min = t
            collapse_idx = i
            collapse_pt = (p1[1] + d1[1] * t, p1[2] + d1[2] * t)
        end
    end
    
    return t_min, collapse_idx, collapse_pt
end

"""
Robust time when p1 + t*d1 == p2 + t*d2 (edge-collapse event).
Returns Inf if (d1-d2) is too small or the fit residual is too large.

This is the 2D system: (p1 - p2) + t*(d1 - d2) = 0
We solve using least-squares and verify consistency with scaled tolerance.
"""
function _ray_ray_intersect_time(p1, d1, p2, d2; tol_rel=1e-10, tol_abs=1e-12)
    # Solve (p1 - p2) + t*(d1 - d2) = 0 in least-squares sense
    ax = d1[1] - d2[1]
    ay = d1[2] - d2[2]
    bx = p2[1] - p1[1]
    by = p2[2] - p1[2]
    
    denom = ax*ax + ay*ay
    denom < tol_abs && return Inf
    
    # Least-squares t
    t = (bx*ax + by*ay) / denom
    t <= 1e-12 && return Inf
    
    # Residual check: p1 + t*d1 should equal p2 + t*d2
    rx = (p1[1] + t*d1[1]) - (p2[1] + t*d2[1])
    ry = (p1[2] + t*d1[2]) - (p2[2] + t*d2[2])
    r = hypot(rx, ry)
    
    # Scale tolerance to problem magnitude (position + displacement scale)
    scale = max(hypot(bx, by), hypot(d1[1]*t, d1[2]*t), hypot(d2[1]*t, d2[2]*t), 1.0)
    tol = max(tol_abs, tol_rel * scale)
    
    return (r <= tol) ? t : Inf
end

"""Build new polygon after batch collapse, also tracking weights."""
function _build_collapsed_polygon_weighted(
    advanced_pts::Vector{NTuple{2,Float64}},
    edge_map::Vector{Int},
    weight_map::Vector{Float64},
    n_active::Int,
    collapse_set::Set{Int},
    collapse_pt_for::Dict{Int, NTuple{2,Float64}}
)
    new_pts = NTuple{2,Float64}[]
    new_edge_map = Int[]
    new_weight_map = Float64[]
    
    for i in 1:n_active
        prev_edge = mod1(i - 1, n_active)
        curr_edge = i
        
        prev_collapsed = prev_edge in collapse_set
        curr_collapsed = curr_edge in collapse_set
        
        if prev_collapsed && curr_collapsed
            continue
        elseif prev_collapsed && !curr_collapsed
            push!(new_pts, collapse_pt_for[prev_edge])
            push!(new_edge_map, edge_map[curr_edge])
            push!(new_weight_map, weight_map[curr_edge])
        elseif !prev_collapsed && curr_collapsed
            continue
        else
            push!(new_pts, advanced_pts[i])
            push!(new_edge_map, edge_map[curr_edge])
            push!(new_weight_map, weight_map[curr_edge])
        end
    end
    
    return new_pts, new_edge_map, new_weight_map
end
