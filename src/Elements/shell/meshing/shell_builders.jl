# ============================================================================
# mesh() - Triangulation only
# ============================================================================

"""
    mesh(corners, n=4; interior_supports=[])

Triangulate a polygon defined by corner nodes. Returns a DelaunayTriangulation object.

# Arguments
- `corners`: Tuple or vector of corner `Node`s in counter-clockwise order (3+ nodes)
- `n`: Subdivision level (default: 4). Higher = finer mesh.
- `interior_supports`: List of interior support lines (Elements or node pairs)

# Example
```julia
tri = mesh((n1, n2, n3, n4), 6)
tri = mesh((n1, n2, n3, n4), 6; interior_supports=[interior_beam])
```
"""
function mesh(
    corners::Union{NTuple{N, Node}, Vector{Node}},
    n::Int = 4;
    interior_supports::Vector = []
) where N
    nc = length(corners)
    @assert nc >= 3 "mesh requires at least 3 corner nodes"
    @assert n >= 1 "Subdivision level must be at least 1"
    
    # Extract corner positions in meters (2D)
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    
    # Convert interior supports to _SupportLine (explicit type avoids Vector{Any} when interior_supports is untyped)
    support_lines = _SupportLine[_to_support_line(s) for s in interior_supports]
    
    # Generate all points (boundary + interior + support lines) without duplicates
    all_points, support_point_indices, _ = _generate_mesh_points_with_supports(boundary_pts, n, support_lines)
    
    # Build constrained triangulation
    boundary_loop = vcat(collect(1:nc), [1])
    return DT.triangulate(all_points; boundary_nodes=[boundary_loop], rng=_mesh_rng())
end

# ============================================================================
# Shell() - Create shell elements
# ============================================================================

"""
    Shell(corners, section; n=4, target_edge_length=nothing, id=:shell, ...)
    Shell(corners, n, section; ...)

Create triangular shell elements from any polygon with support conditions.

# Arguments
- `corners`: Tuple or vector of corner `Node`s in counter-clockwise order
- `section`: `ShellSection` defining thickness and material properties
- `n`: Subdivision level (default: 4). Higher = finer mesh. Ignored when `target_edge_length` is set.

# Keyword Arguments — mesh sizing
- `target_edge_length`: Target element size (default: `0.25u"m"`). Overrides `n` when provided.
  For rectangular panels **without refinement**, uses structured mesh; otherwise Delaunay.
  Set to `nothing` to fall back to legacy `n`-based meshing.
- `mesh_density_warn_shortest_side_m`: Optional. When set (meters), mesh-density warnings use this as the
  reference shortest side instead of the polygon bounding-box `min(Lx, Ly)`. Use for one shell covering
  many structural cells so warnings reflect per-cell `n`, not divisions across the full outline.

# Keyword Arguments — refinement
- `refinement_edge_length`: Element size near refinement targets (e.g., `0.1u"m"`).
  When set, mesh is graded: fine near targets, coarse elsewhere.
- `refinement_radius`: Transition distance from fine to coarse (default: auto = 5× refinement_edge_length).
- `refinement_targets`: Controls where to refine. Default `:auto`.
  - `:auto` — auto-detect from `interior_nodes`, `interior_patches` centroids, and `interior_supports` endpoints.
    Boundary conditions (`edge_support_type`) are NOT considered — mesh refinement and DOF fixity are orthogonal.
  - `:none` — no refinement targets, even if `refinement_edge_length` is set (uniform mesh).
  - `Vector{Node}` — explicit list of points to refine around.

# Keyword Arguments — supports
- `id::Symbol`: Element identifier (default: `:shell`)
- `interior_supports`: List of interior support lines - Elements or (Node, Node) pairs. 
  Creates new mesh nodes with specified fixity along these lines.
- `interior_nodes::Vector{Node}`: Actual Node objects to embed in the mesh. Unlike 
  `interior_supports`, these exact node objects are used (not copies), enabling true 
  structural connectivity with frame elements sharing the same nodes.
- `interior_constrained_lines::Vector{Vector{NTuple{2, Float64}}}`: Interior polylines
  (each a sequence of `(x, y)` points, in meters) whose consecutive-vertex edges become
  **constrained edges** of the Delaunay triangulation. Unlike `interior_supports` (which
  are sampled as points and given boundary conditions), these guarantee triangle edges
  align with the polyline — used for joist / framing lines the slab mesh must respect.
  No boundary conditions are imposed; combine with `interior_nodes` for DOF sharing.
  Empty by default (behavior identical to before when empty).
- `edge_support_type::Symbol`: Fixity for edge nodes (default: `:pinned`)
- `interior_support_type::Symbol`: Fixity for interior support line nodes (default: `:pinned`)

# Support Types
- `:pinned` - x, y, z translations fixed, rotations free (simply-supported)
- `:zfixed` - only z translation fixed (vertical support, free to slide)
- `:fixed` - all 6 DOFs fixed (clamped)
- `:free` - no constraints (for edges that aren't supported)
- `Vector{Bool}` - custom DOF pattern [x, y, z, θx, θy, θz]

# Examples
```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2)

# Target edge length (uniform mesh)
shells = Shell((n1, n2, n3, n4), section; target_edge_length=0.5u"m")

# With graded refinement near columns
shells = Shell((n1, n2, n3, n4), section;
    target_edge_length = 0.5u"m",
    interior_nodes = [col_top],
    refinement_edge_length = 0.1u"m"
)

# Legacy: subdivision level
shells = Shell((n1, n2, n3, n4), 6, section)
```
"""
function Shell(
    corners::Union{NTuple{N, Node}, Vector{Node}},
    n::Int,
    section::ShellSection;
    id::Symbol = :shell,
    interior_supports::Vector = [],
    interior_nodes::Vector{Node} = Node[],
    interior_patches::Vector{ShellPatch} = ShellPatch[],
    interior_constrained_lines::Vector{Vector{NTuple{2, Float64}}} = Vector{NTuple{2, Float64}}[],
    holes::Vector = [],
    edge_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    interior_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    target_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_radius::Union{Quantity, Real, Nothing} = nothing,
    refinement_targets::Union{Symbol, Vector{Node}, Nothing} = :auto,
    mesh_density_warn_shortest_side_m::Union{Float64, Nothing} = nothing,
) where N
    nc = length(corners)
    
    # ===========================================================================
    # Input Validation
    # ===========================================================================
    @assert nc >= 3 "Shell requires at least 3 corner nodes, got $nc"
    
    # Extract corner positions for validation
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    
    hole_loops = [_as_hole_loop(h) for h in holes]
    filter!(h -> !isempty(h), hole_loops)

    area = _shoelace_area(boundary_pts)
    @assert abs(area) > 1e-12 "Shell polygon is degenerate (zero area). Check corner positions."
    @assert !_is_self_intersecting(boundary_pts) "Shell polygon is self-intersecting. Corners must form a simple polygon."
    
    # ===========================================================================
    # Convert parameters to meters
    # ===========================================================================
    h_far = if target_edge_length !== nothing
        target_edge_length isa Quantity ? Float64(ustrip(u"m", target_edge_length)) : Float64(target_edge_length)
    else
        nothing
    end
    
    h_near = if refinement_edge_length !== nothing
        refinement_edge_length isa Quantity ? Float64(ustrip(u"m", refinement_edge_length)) : Float64(refinement_edge_length)
    else
        nothing
    end
    
    r_trans = if refinement_radius !== nothing
        refinement_radius isa Quantity ? Float64(ustrip(u"m", refinement_radius)) : Float64(refinement_radius)
    elseif h_near !== nothing
        5.0 * h_near  # default: transition over 5× refinement size
    else
        nothing
    end
    
    # Convert interior supports (explicit type avoids Vector{Any} when interior_supports is untyped)
    support_lines = _SupportLine[_to_support_line(s) for s in interior_supports]
    
    # ===========================================================================
    # Resolve refinement targets
    # ===========================================================================
    # refinement_targets controls WHERE to refine:
    #   :auto (default) — auto-detect from interior_nodes + patches + support lines
    #   :none           — no refinement targets (uniform mesh even with h_near set)
    #   Vector{Node}    — explicit list of points to refine around
    #   nothing         — treated as :auto (backwards compatibility)
    refine_pts = if h_near !== nothing
        if refinement_targets isa Vector
            # Explicit targets (Vector{Node})
            Tuple{Float64, Float64}[(ustrip(u"m", nd.position[1]), 
                                     ustrip(u"m", nd.position[2])) for nd in refinement_targets]
        elseif refinement_targets === :none
            Tuple{Float64, Float64}[]
        else
            # :auto or nothing — auto-detect from interior features
            _collect_refinement_targets(interior_nodes, support_lines, interior_patches)
        end
    else
        Tuple{Float64, Float64}[]
    end
    
    # Refinement targets along constrained polylines so triangle quality stays
    # acceptable near the constrained edges (only relevant when refinement is on).
    if h_near !== nothing && !isempty(interior_constrained_lines)
        append!(refine_pts, _sample_constrained_line_points(interior_constrained_lines, h_near))
    end
    
    # ===========================================================================
    # Bounding box
    # ===========================================================================
    xmin = minimum(p[1] for p in boundary_pts)
    xmax = maximum(p[1] for p in boundary_pts)
    ymin = minimum(p[2] for p in boundary_pts)
    ymax = maximum(p[2] for p in boundary_pts)
    Lx = xmax - xmin
    Ly = ymax - ymin
    
    # Warn if mesh density looks unreasonable
    if h_far !== nothing
        _warn_mesh_density(h_far, Lx, Ly; reference_shortest_side=mesh_density_warn_shortest_side_m)
    end
    
    # ===========================================================================
    # Decide meshing strategy
    # ===========================================================================
    # Structured mesh only for uniform rectangular panels (no refinement).
    # When refinement is requested, Delaunay + rings is always better because it
    # refines radially/locally instead of creating axis-aligned dense bands.
    # Constrained lines require the Delaunay path (structured T3 meshing cannot
    # honor arbitrary interior constrained edges), so exclude it here.
    use_structured = h_far !== nothing && h_near === nothing && _is_rectangular(boundary_pts) && isempty(support_lines) && isempty(interior_nodes) && isempty(interior_patches) && isempty(interior_constrained_lines) && isempty(hole_loops)
    
    if use_structured
        # ── Structured mesh (T3blockx) ──
        target_xs = [pt[1] for pt in refine_pts]
        target_ys = [pt[2] for pt in refine_pts]
        
        if h_near !== nothing && !isempty(refine_pts)
            xs = _graded_spacing(xmin, xmax, target_xs, h_far, h_near, r_trans)
            ys = _graded_spacing(ymin, ymax, target_ys, h_far, h_near, r_trans)
        else
            xs = _uniform_spacing(xmin, xmax, h_far)
            ys = _uniform_spacing(ymin, ymax, h_far)
        end
        
        return _create_shells_from_structured(
            xs, ys, corners, boundary_pts, interior_nodes, section, id,
            edge_support_type, interior_support_type)
    else
        # ── Delaunay mesh ──
        # Two strategies:
        #   A) Refinement requested (h_near set) → Minimal seed mesh +
        #      Ruppert refinement.  Start with only boundary + interior
        #      nodes + patch vertices.  Let DT.refine! optimally place all
        #      Steiner points — avoids the near-collinear grid points that
        #      cause refine! to crash.
        #   B) No refinement → Dense grid mesh via _generate_mesh_points.
        
        # Compute effective n (used for grid-based path)
        effective_n = if h_far !== nothing
            max(1, ceil(Int, max(Lx, Ly) / h_far))
        else
            n
        end
        @assert effective_n >= 1 "Mesh refinement n must be at least 1"
        
        # Extract interior node positions
        interior_node_pts = Tuple{Float64, Float64}[(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                             for node in interior_nodes]
        
        want_refinement = h_near !== nothing && !isempty(refine_pts)
        
        if want_refinement
            # ── Strategy A: Minimal seed + Ruppert ──
            # Build a minimal point set: boundary corners + interior nodes
            # + patch vertices.  No grid subdivision — refine! handles that.
            tol_dedup = 1e-6
            round_coord(x) = round(Int64, x / tol_dedup)
            
            seed_points = copy(boundary_pts)  # indices 1:nc
            existing_keys = Set{Tuple{Int64,Int64}}(
                (round_coord(p[1]), round_coord(p[2])) for p in boundary_pts)
            
            # Track which seed-point index corresponds to each interior node
            interior_node_map = Dict{Int, Int}()
            for (idx, pt) in enumerate(interior_node_pts)
                key = (round_coord(pt[1]), round_coord(pt[2]))
                if key ∉ existing_keys
                    if _point_inside_polygon(pt, boundary_pts) ||
                       _is_on_boundary_edge(pt, boundary_pts, tol_dedup)
                        push!(seed_points, pt)
                        interior_node_map[idx] = length(seed_points)
                        push!(existing_keys, key)
                    end
                end
            end
            
            # Add patch vertices (no constrained segments — see comment below)
            if !isempty(interior_patches)
                _add_patch_geometry!(seed_points, interior_patches, boundary_pts)
            end
            
            # Add constrained-line vertices and collect their segments. Ruppert
            # refinement may insert Steiner points that subdivide these segments,
            # but the constrained edges themselves are preserved by DT.
            line_segments = isempty(interior_constrained_lines) ? Set{Tuple{Int,Int}}() :
                _add_constrained_line_geometry!(seed_points, interior_constrained_lines, boundary_pts)

            seed_points, seed_remap = _drop_hole_interior(seed_points, nc, hole_loops)
            interior_node_map = Dict(k => seed_remap[v] for (k, v) in interior_node_map if haskey(seed_remap, v))
            line_segments = _remap_pairs(line_segments, seed_remap)
            boundary_nodes = _boundary_curves!(seed_points, nc, hole_loops)
            tri = isempty(line_segments) ?
                DT.triangulate(seed_points; boundary_nodes=boundary_nodes, rng=_mesh_rng()) :
                DT.triangulate(seed_points; boundary_nodes=boundary_nodes, segments=line_segments, rng=_mesh_rng())
            
            # Ruppert refinement: spatially-varying area constraint
            far_h  = h_far !== nothing ? h_far : max(Lx, Ly) / max(effective_n, 1)
            A_far  = 0.433 * far_h^2
            A_near = 0.433 * h_near^2
            r_transition = r_trans !== nothing ? r_trans : 5.0 * h_near

            refine_grid = _PointGrid(refine_pts)
            function _needs_refine(_tri, T)
                i, j, k = DT.triangle_vertices(T)
                p, q, r = DT.get_point(_tri, i, j, k)
                A = abs(DT.triangle_area(p, q, r))
                cx = (p[1] + q[1] + r[1]) / 3
                cy = (p[2] + q[2] + r[2]) / 3
                d_min = _min_dist(refine_grid, cx, cy)
                α = clamp(d_min / r_transition, 0.0, 1.0)
                A_limit = A_near + α * (A_far - A_near)
                return A > A_limit
            end

            max_pts = max(5000, 20 * length(seed_points))
            tri_backup = deepcopy(tri)
            try
                DT.refine!(tri; min_angle=30.0, custom_constraint=_needs_refine,
                           max_points=max_pts, rng=_mesh_rng())
            catch e
                @warn "Ruppert refinement failed — falling back to grid mesh" exception=(e, catch_backtrace())
                # Fall back to grid-based mesh (Strategy B)
                tri = nothing  # signal to use grid path below
            end
            
            if tri !== nothing
                # Build node map from the refined triangulation
                z_val = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
                points = DT.get_points(tri)
                node_map = Dict{Int, Node}()
                edge_node_indices = Set{Int}()
                
                tol = 1e-6
                for (idx, pt) in enumerate(points)
                    if _is_on_boundary_edge(pt, boundary_pts, tol)
                        push!(edge_node_indices, idx)
                    end
                end
                
                for (i, corner) in enumerate(corners)
                    node_map[i] = corner
                end
                
                for (input_idx, point_idx) in interior_node_map
                    node_map[point_idx] = interior_nodes[input_idx]
                end
                
                for i in (nc+1):length(points)
                    haskey(node_map, i) && continue
                    pt = points[i]
                    pos = [pt[1], pt[2], z_val] .* u"m"
                    fixity = if _on_hole_edge(pt, hole_loops, tol)
                        :free
                    elseif i in edge_node_indices
                        edge_support_type
                    else
                        :free
                    end
                    node_map[i] = _create_node_with_fixity(pos, fixity, id)
                end
                
                shells = _create_shells_from_tri(tri, node_map, section, id)
                
                if !isempty(interior_patches)
                    _apply_patches!(shells, interior_patches)
                end
                
                return shells
            end
            # If we fall through here, tri was set to nothing — use grid path
        end
        
        # ── Strategy B: Dense grid mesh (no Ruppert) ──
        # Pass constrained lines so grid/edge/support points landing on a
        # constrained segment interior are dropped — otherwise DT would split
        # the constrained edge at those collinear points.
        all_points, support_point_indices, interior_node_map = _generate_mesh_points_with_supports(
            boundary_pts, effective_n, support_lines, interior_node_pts, interior_constrained_lines)
        
        # Add patch vertices as regular mesh points (no constrained segments).
        # _apply_patches! assigns stiffened sections by centroid check, so
        # strict edge conformity is not required.
        if !isempty(interior_patches)
            _add_patch_geometry!(all_points, interior_patches, boundary_pts)
        end
        
        # Add constrained-line vertices and collect their segments for a
        # constrained Delaunay triangulation. With no interior collinear points
        # (dropped above), each consecutive-vertex pair survives as a mesh edge.
        line_segments = isempty(interior_constrained_lines) ? Set{Tuple{Int,Int}}() :
            _add_constrained_line_geometry!(all_points, interior_constrained_lines, boundary_pts)

        all_points, grid_remap = _drop_hole_interior(all_points, nc, hole_loops)
        support_point_indices = Set(grid_remap[i] for i in support_point_indices if haskey(grid_remap, i))
        interior_node_map = Dict(k => grid_remap[v] for (k, v) in interior_node_map if haskey(grid_remap, v))
        line_segments = _remap_pairs(line_segments, grid_remap)
        boundary_nodes = _boundary_curves!(all_points, nc, hole_loops)
        tri = isempty(line_segments) ?
            DT.triangulate(all_points; boundary_nodes=boundary_nodes, rng=_mesh_rng()) :
            DT.triangulate(all_points; boundary_nodes=boundary_nodes, segments=line_segments, rng=_mesh_rng())
        
        z_val = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
        
        points = DT.get_points(tri)
        node_map = Dict{Int, Node}()
        edge_node_indices = Set{Int}()
        
        tol = 1e-6
        for (idx, pt) in enumerate(points)
            if _is_on_boundary_edge(pt, boundary_pts, tol)
                push!(edge_node_indices, idx)
            end
        end
        
        for (i, corner) in enumerate(corners)
            node_map[i] = corner
        end
        
        for (input_idx, point_idx) in interior_node_map
            node_map[point_idx] = interior_nodes[input_idx]
        end
        
        for i in (nc+1):length(points)
            haskey(node_map, i) && continue
            pt = points[i]
            pos = [pt[1], pt[2], z_val] .* u"m"
            
            if _on_hole_edge(pt, hole_loops, tol)
                fixity = :free
            elseif i in support_point_indices
                fixity = interior_support_type
            elseif i in edge_node_indices
                fixity = edge_support_type
            else
                fixity = :free
            end
            
            node_map[i] = _create_node_with_fixity(pos, fixity, id)
        end
        
        shells = _create_shells_from_tri(tri, node_map, section, id)
        
        if !isempty(interior_patches)
            _apply_patches!(shells, interior_patches)
        end
        
        return shells
    end
end

# Keyword-n version (default target_edge_length=0.25m)
function Shell(
    corners::Union{NTuple{N, Node}, Vector{Node}},
    section::ShellSection;
    n::Int = 4,
    id::Symbol = :shell,
    interior_supports::Vector = [],
    interior_nodes::Vector{Node} = Node[],
    interior_patches::Vector{ShellPatch} = ShellPatch[],
    interior_constrained_lines::Vector{Vector{NTuple{2, Float64}}} = Vector{NTuple{2, Float64}}[],
    holes::Vector = [],
    edge_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    interior_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    target_edge_length::Union{Quantity, Real, Nothing} = 0.25u"m",
    refinement_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_radius::Union{Quantity, Real, Nothing} = nothing,
    refinement_targets::Union{Symbol, Vector{Node}, Nothing} = :auto,
    mesh_density_warn_shortest_side_m::Union{Float64, Nothing} = nothing,
) where N
    return Shell(corners, n, section; 
        id=id, 
        interior_supports=interior_supports,
        interior_nodes=interior_nodes,
        interior_patches=interior_patches,
        interior_constrained_lines=interior_constrained_lines,
        holes=holes,
        edge_support_type=edge_support_type,
        interior_support_type=interior_support_type,
        target_edge_length=target_edge_length,
        refinement_edge_length=refinement_edge_length,
        refinement_radius=refinement_radius,
        refinement_targets=refinement_targets,
        mesh_density_warn_shortest_side_m=mesh_density_warn_shortest_side_m,
    )
end
