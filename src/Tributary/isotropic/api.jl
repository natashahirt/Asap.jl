# =============================================================================
# Straight Skeleton using DCEL (Weighted Tributary Areas)
# =============================================================================

"""
    get_tributary_polygons_isotropic(vertices; weights=nothing, holes=nothing)

Compute tributary polygons using straight skeleton with DCEL data structure.

## Arguments
- `vertices::Vector{<:Point}`: Polygon vertices as Meshes.Point objects (any Unitful length - 
  automatically converted to meters internally)
- `weights::Union{Nothing, AbstractVector{<:Real}}`: Optional edge weights (one per edge).
  Higher weight = faster shrink = smaller tributary area. Default `nothing` = all weights 1.0.
- `holes`: Optional openings, each a vertex loop (`Point`s or `(x, y)` in meters).
  An opening edge is not a support. The skeleton grows from the supported
  edges only, and each edge's region loses the part inside an opening.

## Returns
`Vector{TributaryPolygon}` in parametric form. All length values (`d`, `area`) are in meters.
Use `vertices(trib, beam_start, beam_end)` with beam coords in meters to get absolute coords.
Areas are net of openings and sum to the net floor area; this throws otherwise.

At a re-entrant corner an edge's region runs past the edge's end. That part
is folded back onto the edge about the end ([`_edge_slice_profile`](@ref)).
The wavefront has edge events only, with no split events, so the regions of a
concave face can differ from the true straight skeleton. They still tile the
face, and the balance check still holds.
"""
function get_tributary_polygons_isotropic(
    vertices::Vector{<:Point};
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    holes = nothing
)
    m = length(vertices)  # Original vertex count
    m >= 3 || return TributaryPolygon[]
    
    # Convert to simple 2D coords in METERS
    pts_orig = [_to_2d(v) for v in vertices]
    
    # Ensure CCW orientation (algorithm assumes interior is on LEFT)
    pts_orig = _ensure_ccw(pts_orig)
    original_pts = copy(pts_orig)  # Keep for area calculation
    _, hole_loops = _span_region(pts_orig, holes)
    net_area = _net_area(original_pts, hole_loops)
    
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
    
    # Original edges of each simplified edge k, in chain order. A chain that
    # wraps past vertex 1 lists its last edges first, which is the order of
    # `split_collinear_face`'s pieces.
    groups = [Int[] for _ in 1:n]  # n = simplified edge count
    for k in 1:n
        i = keep_idx[k]
        b = keep_idx[mod1(k + 1, n)]
        while i != b
            push!(groups[k], i)
            i = mod1(i + 1, m)
        end
    end
    
    # Map results back to original edges, using geometric splitting for collinear chains
    results = TributaryPolygon[]
    piece(i, poly, area) = _make_tributary_with_depths(i, poly, pts_orig, area, 0.0; holes = hole_loops)
    
    # Process each simplified edge k and its corresponding original edges
    for k in 1:n
        idxs = groups[k]
        isempty(idxs) && continue
        
        poly, area_k = simp_results[k]
        if isempty(poly) || area_k <= 0
            # No valid polygon for this simplified edge
            for i in idxs
                push!(results, piece(i, NTuple{2,Float64}[], 0.0))
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
            push!(results, piece(idxs[1], poly, area_k))
        else
            # One piece per sub-edge, cut at the chain's interior vertices.
            split_polys = split_collinear_face(poly, chain_pts)
            for (j, i) in enumerate(idxs)
                split_poly = split_polys[j]
                push!(results, piece(i, split_poly, abs(_polygon_area(split_poly))))
            end
        end
    end
    
    # Handle any original edges that didn't map to a simplified edge
    # (should be rare, but handle for completeness)
    processed_edges = Set(r.local_edge_idx for r in results)
    for i in 1:m
        if i ∉ processed_edges && (orig_to_simp[i] == 0 || orig_to_simp[i] > n)
            push!(results, piece(i, NTuple{2,Float64}[], 0.0))
        end
    end
    
    # Sort results by edge index
    sort!(results, by=r -> r.local_edge_idx)
    results = [_with_fraction(r, net_area > 0 ? r.area / net_area : 0.0) for r in results]
    assert_tributary_balance(results, pts_orig, net_area; context = " (two-way)")
    return results
end

