# =============================================================================
# Straight Skeleton using DCEL (Weighted Tributary Areas)
# =============================================================================

"""
    get_tributary_polygons_isotropic(vertices; weights=nothing)

Compute tributary polygons using straight skeleton with DCEL data structure.

## Arguments
- `vertices::Vector{<:Point}`: Polygon vertices as Meshes.Point objects (any Unitful length - 
  automatically converted to meters internally)
- `weights::Union{Nothing, AbstractVector{<:Real}}`: Optional edge weights (one per edge).
  Higher weight = faster shrink = smaller tributary area. Default `nothing` = all weights 1.0.

## Returns
`Vector{TributaryPolygon}` in parametric form. All length values (`d`, `area`) are in meters.
Use `vertices(trib, beam_start, beam_end)` with beam coords in meters to get absolute coords.
"""
function get_tributary_polygons_isotropic(
    vertices::Vector{<:Point};
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing
)
    m = length(vertices)  # Original vertex count
    m >= 3 || return TributaryPolygon[]
    
    # Convert to simple 2D coords in METERS
    pts_orig = [_to_2d(v) for v in vertices]
    
    # Ensure CCW orientation (algorithm assumes interior is on LEFT)
    pts_orig = _ensure_ccw(pts_orig)
    original_pts = copy(pts_orig)  # Keep for area calculation
    
    # Handle weights
    if isnothing(weights)
        weights_orig = ones(m)
    else
        length(weights) == m || error("weights must have same length as vertices ($m)")
        weights_orig = Float64.(weights)
        all(w -> w > 0, weights_orig) || error("all weights must be positive")
    end
    
    # Simplify collinear vertices (removes 180° degeneracies)
    pts, keep_idx = simplify_collinear_polygon(pts_orig; tol=1e-12)
    n = length(pts)
    n >= 3 || return TributaryPolygon[]  # Degenerate after simplification
    
    # Build mapping: orig_to_simp[i] = simplified edge that contains original edge i
    # Also compute simplified weights (average weight of merged edges)
    orig_to_simp = fill(0, m)
    n_s = length(keep_idx)
    simp_weights = zeros(n_s)
    simp_edge_counts = zeros(Int, n_s)
    
    for k in 1:n_s
        a = keep_idx[k]
        b = keep_idx[mod1(k + 1, n_s)]
        i = a
        while i != b
            orig_to_simp[i] = k
            simp_weights[k] += weights_orig[i]
            simp_edge_counts[k] += 1
            i = mod1(i + 1, m)
        end
    end
    
    # Average the weights for merged edges
    for k in 1:n_s
        if simp_edge_counts[k] > 0
            simp_weights[k] /= simp_edge_counts[k]
        else
            simp_weights[k] = 1.0
        end
    end
    
    # Initialize DCEL and vertex registry
    dcel = DCEL()
    registry = VertexRegistry(dcel; tol=1e-9)
    
    # =========================================================================
    # Step 1: Build boundary halfedges for the simplified polygon
    # =========================================================================
    boundary_vertices = Int[]
    for i in 1:n
        v_idx = get_or_create_vertex!(registry, pts[i])
        push!(boundary_vertices, v_idx)
    end
    
    # Create boundary halfedge pairs
    # Edge i (from vertex i to vertex i+1) has:
    #   - Inner halfedge with face = i (tributary face for simplified edge i)
    #   - Outer halfedge with face = 0 (outside)
    boundary_inner = Int[]  # Inner halfedges (one per simplified edge)
    boundary_outer = Int[]  # Outer halfedges
    
    for i in 1:n
        v1 = boundary_vertices[i]
        v2 = boundary_vertices[mod1(i + 1, n)]
        h_inner, h_outer = create_halfedge_pair!(dcel, v1, v2, i, 0)
        push!(boundary_inner, h_inner)
        push!(boundary_outer, h_outer)
    end
    
    # =========================================================================
    # Step 2: Wavefront propagation with skeleton arc recording
    # =========================================================================
    
    # Current wavefront state
    current_pts = copy(pts)
    n_active = n
    
    # edge_map[i] = which original edge the current wavefront edge i represents
    edge_map = collect(1:n)
    
    # weight_map[i] = weight of current wavefront edge i
    weight_map = copy(simp_weights)
    
    # Check if isotropic (all weights equal to each other)
    # For isotropic case, all weights should be 1.0, but we check if they're all equal
    is_isotropic = length(weight_map) < 2 || all(w -> abs(w - weight_map[1]) < 1e-10, weight_map)
    
    # Track which original edges are still active
    edge_active = trues(n)
    
    while n_active > 2
        if is_isotropic
            # Isotropic case: use angle bisector approach
            bisectors, speeds = _compute_bisectors(current_pts, n_active)
            
            # Find all edges collapsing at the next event time
            t_min, collapses = _find_all_collapses(current_pts, n_active, bisectors, speeds)
            
            if t_min == Inf || t_min <= 1e-10 || isempty(collapses)
                break
            end
            
            # Advance all vertices to time t_min using bisector approach
            new_pts = Vector{NTuple{2,Float64}}(undef, n_active)
            for i in 1:n_active
                bx, by = bisectors[i]
                s = speeds[i]
                px, py = current_pts[i]
                new_pts[i] = (px + bx * s * t_min, py + by * s * t_min)
            end
        else
            # Weighted case: use affine motion model
            # x0 is current_pts (t=0 state), we only compute velocities v
            normals, _ = _edge_lines_ccw(current_pts, weight_map)
            v = _vertex_affine_velocity(normals, weight_map)
            
            # Filter out degenerate vertices (NaN velocities)
            valid = [i for i in 1:n_active if isfinite(v[i][1]) && isfinite(v[i][2])]
            if length(valid) < 3
                break  # Too few valid vertices
            end
            
            # Find all edges collapsing at the next event time
            x0 = current_pts  # Use current positions as t=0 state
            t_min, collapses = _find_all_collapses_affine(x0, v, n_active)
            
            if t_min == Inf || t_min <= 1e-10 || isempty(collapses)
                break
            end
            
            # Advance all vertices to time t_min using affine motion: x(t) = x0 + v*t
            new_pts = Vector{NTuple{2,Float64}}(undef, n_active)
            for i in 1:n_active
                if isfinite(v[i][1]) && isfinite(v[i][2])
                    new_pts[i] = (x0[i][1] + v[i][1]*t_min, x0[i][2] + v[i][2]*t_min)
                else
                    # Degenerate vertex - keep current position
                    new_pts[i] = current_pts[i]
                end
            end
        end
        
        # Record skeleton arcs for each active vertex trajectory
        # Each vertex i sits between wavefront edges (i-1) and i
        # So its trajectory separates faces edge_map[mod1(i-1, n_active)] and edge_map[i]
        for i in 1:n_active
            p_old = current_pts[i]
            p_new = new_pts[i]
            
            # Skip degenerate arcs
            if _dist(p_old, p_new) < 1e-10
                continue
            end
            
            # Adjacent faces (original edge IDs)
            face_left = edge_map[mod1(i - 1, n_active)]
            face_right = edge_map[i]
            
            # Record the skeleton arc
            _record_skeleton_arc!(dcel, registry, p_old, p_new, face_left, face_right)
        end
        
        # Build lookup structures for batch processing
        collapse_set = Set(c[1] for c in collapses)
        collapse_pt_for = Dict(c[1] => c[2] for c in collapses)
        
        # Mark collapsing edges
        for (idx, _) in collapses
            edge_active[edge_map[idx]] = false
        end
        
        # Handle case where ALL edges collapse (final convergence)
        if length(collapse_set) == n_active
            # All edges collapse simultaneously
            cps = [pt for (_, pt) in collapses]
            
            # Check if all collapse points are coincident (within tolerance)
            # For isotropic symmetric cases they coincide; with weights they may form a segment
            tol_meet = 1e-8
            all_coincident = true
            for i in 2:length(cps)
                if _dist(cps[1], cps[i]) > tol_meet
                    all_coincident = false
                    break
                end
            end
            
            if all_coincident
                # All points coincide — use average (or just first point)
                meet = (sum(p[1] for p in cps) / length(cps), 
                        sum(p[2] for p in cps) / length(cps))
                
                # Record arcs from each vertex to meeting point
                for i in 1:n_active
                    p_old = current_pts[i]
                    if _dist(p_old, meet) > 1e-10
                        face_left = edge_map[mod1(i - 1, n_active)]
                        face_right = edge_map[i]
                        _record_skeleton_arc!(dcel, registry, p_old, meet, face_left, face_right)
                    end
                end
            else
                # Collapse points form a segment (roof ridge) — record to segment endpoints
                # Find bounding box of collapse points
                xs = [p[1] for p in cps]
                ys = [p[2] for p in cps]
                x_min, x_max = extrema(xs)
                y_min, y_max = extrema(ys)
                
                # Use endpoints of the bounding box (or cluster and use cluster centers)
                # For now, use min/max points as segment endpoints
                meet1 = (x_min, y_min)
                meet2 = (x_max, y_max)
                
                # Record arcs to both endpoints (or to nearest endpoint)
                for i in 1:n_active
                    p_old = current_pts[i]
                    face_left = edge_map[mod1(i - 1, n_active)]
                    face_right = edge_map[i]
                    
                    # Choose closer endpoint
                    d1 = _dist(p_old, meet1)
                    d2 = _dist(p_old, meet2)
                    meet = d1 < d2 ? meet1 : meet2
                    
                    if _dist(p_old, meet) > 1e-10
                        _record_skeleton_arc!(dcel, registry, p_old, meet, face_left, face_right)
                    end
                end
            end
            break
        end
        
        # Build new polygon handling batch collapses (also updates weights)
        new_current_pts, new_edge_map, new_weight_map = _build_collapsed_polygon_weighted(
            new_pts, edge_map, weight_map, n_active, collapse_set, collapse_pt_for
        )
        
        # Clean up degenerate edges/vertices introduced by collapse rebuild
        current_pts, edge_map, weight_map = cleanup_wavefront(
            new_current_pts, new_edge_map, new_weight_map;
            eps_len=1e-9, eps_col=1e-12
        )
        n_active = length(current_pts)
        
        # Recheck if still isotropic (after collapses, weights might have changed)
        is_isotropic = length(weight_map) < 2 || all(w -> abs(w - weight_map[1]) < 1e-10, weight_map)
    end
    
    # =========================================================================
    # Step 3: Final convergence - remaining vertices meet at center
    # =========================================================================
    if n_active >= 2
        # Check if still isotropic (all weights equal)
        is_iso_final = length(weight_map) < 2 || all(w -> abs(w - weight_map[1]) < 1e-10, weight_map)
        
        if is_iso_final
            # Isotropic case
            bisectors, speeds = _compute_bisectors(current_pts, n_active)
            t_min, _, final_pt = _find_next_collapse_single(current_pts, n_active, bisectors, speeds)
        else
            # Weighted case: use affine motion
            normals, _ = _edge_lines_ccw(current_pts, weight_map)
            v = _vertex_affine_velocity(normals, weight_map)
            x0 = current_pts  # Use current positions as t=0 state
            
            # Find next collapse
            t_min = Inf
            final_pt = (0.0, 0.0)
            for i in 1:n_active
                next_i = mod1(i + 1, n_active)
                if isfinite(v[i][1]) && isfinite(v[i][2]) && isfinite(v[next_i][1]) && isfinite(v[next_i][2])
                    t = _meet_time_affine(x0[i], v[i], x0[next_i], v[next_i])
                    if isfinite(t) && t > 1e-10 && t < t_min
                        t_min = t
                        final_pt = (x0[i][1] + v[i][1]*t, x0[i][2] + v[i][2]*t)
                    end
                end
            end
        end
        
        if t_min < Inf && t_min > 1e-10
            # Record final skeleton arcs to the center point
            for i in 1:n_active
                p_old = current_pts[i]
                if _dist(p_old, final_pt) > 1e-10
                    face_left = edge_map[mod1(i - 1, n_active)]
                    face_right = edge_map[i]
                    _record_skeleton_arc!(dcel, registry, p_old, final_pt, face_left, face_right)
                end
            end
        end
    end
    
    # =========================================================================
    # Step 4: Compute next/prev pointers using DCEL rotation system
    # =========================================================================
    compute_next_prev!(dcel)
    
    # =========================================================================
    # Step 4b: Insert artificial bisectors (micro-spokes, not self-loops)
    # =========================================================================
    n_bisectors = insert_artificial_bisectors!(dcel, registry)
    if n_bisectors > 0
        # Recompute next/prev after inserting bisectors
        compute_next_prev!(dcel)
    end
    
    # =========================================================================
    # Step 4c: Validate DCEL
    # =========================================================================
    valid, errs = validate_dcel(dcel)
    valid || @warn(join(errs, "\n"))
    
    # =========================================================================
    # Step 5: Extract tributary polygons for simplified faces
    # =========================================================================
    total_area = abs(_polygon_area(original_pts))
    simp_results = Vector{Tuple{Vector{NTuple{2,Float64}}, Float64}}(undef, n)
    
    for i in 1:n
        # Use cycle walking as primary (rotation system ensures correctness)
        poly_verts = extract_face_polygon(dcel, i)
        
        if isempty(poly_verts) || length(poly_verts) < 3
            # Fallback: try edge-based extraction for debugging
            poly_verts = extract_face_polygon_by_edges(dcel, i)
        end
        
        if isempty(poly_verts) || length(poly_verts) < 3
            simp_results[i] = (NTuple{2,Float64}[], 0.0)
        else
            area = abs(_polygon_area(poly_verts))
            simp_results[i] = (poly_verts, area)
        end
    end
    
    # =========================================================================
    # Step 6: Map simplified results back to original edges
    # =========================================================================
    
    # Build groups of original edges that map to each simplified edge k
    groups = [Int[] for _ in 1:n]  # n = simplified edge count
    for i in 1:m
        k = orig_to_simp[i]
        if 1 <= k <= n
            push!(groups[k], i)
        end
    end
    
    # Compute original edge lengths (in original polygon indexing)
    function _edge_len(orig_pts, i)
        p = orig_pts[i]
        q = orig_pts[mod1(i+1, length(orig_pts))]
        return hypot(q[1]-p[1], q[2]-p[2])
    end
    orig_edge_len = [_edge_len(pts_orig, i) for i in 1:m]
    
    # Map results back to original edges, using geometric splitting for collinear chains
    results = TributaryPolygon[]
    
    # Process each simplified edge k and its corresponding original edges
    for k in 1:n
        idxs = groups[k]
        isempty(idxs) && continue
        
        poly, area_k = simp_results[k]
        if isempty(poly) || area_k <= 0
            # No valid polygon for this simplified edge
            for i in idxs
                push!(results, _make_tributary_with_depths(i, NTuple{2,Float64}[], pts_orig, 0.0, 0.0))
            end
            continue
        end
        
        # Extract chain of vertices from original polygon for this collinear chain
        # The simplified edge k corresponds to original vertices from keep_idx[k] to keep_idx[k+1]
        a = keep_idx[k]
        b = keep_idx[mod1(k + 1, n)]
        
        # Build chain of vertices (in original polygon order)
        chain_pts = NTuple{2,Float64}[]
        i_chain = a
        while true
            push!(chain_pts, pts_orig[i_chain])
            if i_chain == b
                break
            end
            i_chain = mod1(i_chain + 1, m)
        end
        
        if length(idxs) == 1
            # Single edge: no splitting needed, use full polygon
            i = idxs[1]
            frac_i = total_area > 0 ? area_k / total_area : 0.0
            push!(results, _make_tributary_with_depths(i, poly, pts_orig, area_k, frac_i))
        else
            # Multiple edges: split polygon geometrically
            split_polys = split_collinear_face(poly, chain_pts)
            
            # Map each split polygon to its corresponding original edge
            for (j, i) in enumerate(idxs)
                if j <= length(split_polys)
                    split_poly = split_polys[j]
                    area_i = abs(_polygon_area(split_poly))
                    frac_i = total_area > 0 ? area_i / total_area : 0.0
                    push!(results, _make_tributary_with_depths(i, split_poly, pts_orig, area_i, frac_i))
                else
                    # Fallback: proportional area split if splitting failed
                    denom = sum(orig_edge_len[idx] for idx in idxs)
                    share = denom > 0 ? orig_edge_len[i] / denom : 1/length(idxs)
                    area_i = area_k * share
                    frac_i = total_area > 0 ? area_i / total_area : 0.0
                    push!(results, _make_tributary_with_depths(i, poly, pts_orig, area_i, frac_i))
                end
            end
        end
    end
    
    # Handle any original edges that didn't map to a simplified edge
    # (should be rare, but handle for completeness)
    processed_edges = Set(r.local_edge_idx for r in results)
    for i in 1:m
        if i ∉ processed_edges && (orig_to_simp[i] == 0 || orig_to_simp[i] > n)
            push!(results, _make_tributary_with_depths(i, NTuple{2,Float64}[], pts_orig, 0.0, 0.0))
        end
    end
    
    # Sort results by edge index
    sort!(results, by=r -> r.local_edge_idx)
    
    return results
end

# =============================================================================
# Helper Functions (DCEL-specific versions)
# =============================================================================

"""Distance between two points."""
_dist(p1, p2) = hypot(p1[1] - p2[1], p1[2] - p2[2])

"""
Clip polygon poly by half-plane: dot(x - p0, n) >= 0
Sutherland–Hodgman clipping against a half-plane.
"""
function clip_halfplane(poly::Vector{NTuple{2,Float64}},
                        p0::NTuple{2,Float64},
                        n::NTuple{2,Float64};
                        eps=1e-12)
    isempty(poly) && return poly

    function inside(p)
        ((p[1]-p0[1])*n[1] + (p[2]-p0[2])*n[2]) >= -eps
    end

    function intersect(a, b)
        # find t in [0,1] where dot((a + t(b-a)) - p0, n) == 0
        ax = (a[1]-p0[1])*n[1] + (a[2]-p0[2])*n[2]
        bx = (b[1]-p0[1])*n[1] + (b[2]-p0[2])*n[2]
        denom = bx - ax
        if abs(denom) < eps
            return b  # nearly parallel; fallback
        end
        t = -ax / denom
        t = clamp(t, 0.0, 1.0)
        return (a[1] + t*(b[1]-a[1]), a[2] + t*(b[2]-a[2]))
    end

    out = NTuple{2,Float64}[]
    nV = length(poly)
    for i in 1:nV
        a = poly[i]
        b = poly[mod1(i+1, nV)]
        ina = inside(a)
        inb = inside(b)
        if ina && inb
            push!(out, b)
        elseif ina && !inb
            push!(out, intersect(a,b))
        elseif !ina && inb
            push!(out, intersect(a,b))
            push!(out, b)
        end
    end

    return out
end

"""
Split a merged-face polygon P for a collinear chain p0..pk into per-subedge polygons.
Returns a vector of polygons, one for each sub-edge in the chain.
"""
function split_collinear_face(P::Vector{NTuple{2,Float64}},
                              chain_pts::Vector{NTuple{2,Float64}})
    k = length(chain_pts) - 1
    k >= 1 || return [P]

    # tangent along the chain
    tvec = (chain_pts[end][1]-chain_pts[1][1], chain_pts[end][2]-chain_pts[1][2])
    tl = hypot(tvec...)
    tl < 1e-15 && error("degenerate chain")
    t = (tvec[1]/tl, tvec[2]/tl)

    polys = Vector{Vector{NTuple{2,Float64}}}(undef, k)
    for i in 1:k
        Pi = P
        # keep points with dot(x - p_i, t) >= 0
        Pi = clip_halfplane(Pi, chain_pts[i], t)
        # keep points with dot(x - p_{i+1}, t) <= 0  <=> dot(x - p_{i+1}, -t) >= 0
        Pi = clip_halfplane(Pi, chain_pts[i+1], (-t[1], -t[2]))
        polys[i] = Pi
    end
    return polys
end

"""Record a skeleton arc from p_old to p_new separating face_left and face_right."""
function _record_skeleton_arc!(dcel::DCEL, registry::VertexRegistry, 
                                p_old::NTuple{2,Float64}, p_new::NTuple{2,Float64},
                                face_left::Int, face_right::Int)
    v1 = get_or_create_vertex!(registry, p_old)
    v2 = get_or_create_vertex!(registry, p_new)
    
    # Don't create self-loops
    v1 == v2 && return
    
    # Create halfedge pair marked as SKELETON edges (face separators):
    # Walking from v1 to v2, face_left is on the LEFT
    # DCEL convention: halfedge's face = the face on its LEFT
    # h1: v1 → v2 with face = face_left
    # h2: v2 → v1 with face = face_right (original right is now on left)
    # These edges SEPARATE faces and should not be crossed during face walks
    create_halfedge_pair!(dcel, v1, v2, face_left, face_right; is_skeleton=true)
end

"""Compute bisectors and speeds for active vertices (isotropic, all weights = 1)."""
function _compute_bisectors(pts::Vector{NTuple{2,Float64}}, n::Int)
    bisectors = Vector{NTuple{2,Float64}}(undef, n)
    speeds = Vector{Float64}(undef, n)
    
    for i in 1:n
        prev_i = mod1(i - 1, n)
        next_i = mod1(i + 1, n)
        
        p_prev = pts[prev_i]
        p_curr = pts[i]
        p_next = pts[next_i]
        
        v_in = (p_curr[1] - p_prev[1], p_curr[2] - p_prev[2])
        v_out = (p_next[1] - p_curr[1], p_next[2] - p_curr[2])
        
        len_in = hypot(v_in...)
        len_out = hypot(v_out...)
        
        if len_in < 1e-10 || len_out < 1e-10
            bisectors[i] = (0.0, 0.0)
            speeds[i] = 0.0
            continue
        end
        
        v_in = v_in ./ len_in
        v_out = v_out ./ len_out
        
        # Inward normals (90° CCW rotation)
        n_in = (-v_in[2], v_in[1])
        n_out = (-v_out[2], v_out[1])
        
        bx, by = n_in[1] + n_out[1], n_in[2] + n_out[2]
        blen = hypot(bx, by)
        
        if blen < 1e-10
            bisectors[i] = (0.0, 0.0)
            speeds[i] = 0.0
        else
            bisectors[i] = (bx / blen, by / blen)
            speeds[i] = 2.0 / blen
        end
    end
    
    return bisectors, speeds
end

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
