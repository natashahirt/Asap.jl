# VaultShell - Curved parabolic vault shell mesh
# ============================================================================

"""
    VaultShell(corners, section, span_axis, rise; kwargs...)

Create triangular shell elements on a curved parabolic vault surface.

The vault has a parabolic cross-section along the span direction:
    z(s) = z₀ + (4h/L²) · s · (L - s)

where L is the span length and h is the rise at crown.

Uses 2D Delaunay triangulation (with optional Ruppert refinement) projected onto
the parabolic surface. This gives the same refinement behavior as flat `Shell`
near columns/supports.

# Arguments
- `corners::NTuple{4, Node}`: Four corner nodes defining the rectangular plan (CCW order)
- `section::ShellSection`: Shell thickness and material properties
- `span_axis::Tuple{Float64, Float64}`: Unit vector along span direction (vault arches in this direction)
- `rise::Real`: Rise at crown (in meters or Unitful Length)

# Keyword Arguments — mesh sizing
- `target_edge_length`: Target element size (default: `0.25u"m"`). Overrides `n` when provided.
- `n::Int=20`: Subdivision level (used when target_edge_length is nothing)

# Keyword Arguments — refinement
- `interior_nodes::Vector{Node}`: Nodes to embed in mesh (e.g., column tops for structural connectivity)
- `refinement_edge_length`: Element size near refinement targets (e.g., `0.1u"m"`)
- `refinement_radius`: Transition distance from fine to coarse (default: 5× refinement_edge_length)
- `refinement_targets`: `:auto` (from interior_nodes), `:none`, or `Vector{Node}`
- `mesh_density_warn_shortest_side_m`: Optional meters; same as flat `Shell` — use for per-cell reference when the plan outline is much larger than a typical cell.

# Keyword Arguments — supports
- `id::Symbol=:vault`: Element identifier
- `edge_support_type::Symbol=:free`: Fixity for edge nodes

# Returns
`Vector{ShellTri3}`: Flat-facet triangular elements approximating the curved vault surface

# Example
```julia
n1 = Node([0.0u"m", 0.0u"m", 3.0u"m"], :pinned)
n2 = Node([6.0u"m", 0.0u"m", 3.0u"m"], :pinned)
n3 = Node([6.0u"m", 4.0u"m", 3.0u"m"], :pinned)
n4 = Node([0.0u"m", 4.0u"m", 3.0u"m"], :pinned)

section = ShellSection(0.10u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")

# Basic vault
shells = VaultShell((n1, n2, n3, n4), section, (1.0, 0.0), 0.75u"m")

# With refinement near column
shells = VaultShell((n1, n2, n3, n4), section, (1.0, 0.0), 0.75u"m";
    interior_nodes = [col_top],
    refinement_edge_length = 0.1u"m")
```
"""
function VaultShell(
    corner_input::Union{Tuple, AbstractVector},
    section::ShellSection,
    span_axis::Tuple{<:Real, <:Real},
    rise::Union{Quantity, Real};
    n::Int = 20,
    id::Symbol = :vault,
    interior_nodes::Vector{Node} = Node[],
    edge_support_type::Union{Symbol, Vector{Bool}} = :free,
    target_edge_length::Union{Quantity, Real, Nothing} = 0.25u"m",
    refinement_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_radius::Union{Quantity, Real, Nothing} = nothing,
    refinement_targets::Union{Symbol, Vector{Node}, Nothing} = :auto,
    mesh_density_warn_shortest_side_m::Union{Float64, Nothing} = nothing,
)
    length(corner_input) == 4 || throw(ArgumentError(
        "VaultShell requires exactly 4 corner nodes; received $(length(corner_input))."))
    all(node -> node isa Node, corner_input) || throw(ArgumentError(
        "VaultShell corners must all be Node objects."))
    corners = Node[node for node in corner_input]
    nc = length(corners)
    
    # Convert rise to meters
    rise_m = rise isa Quantity ? Float64(ustrip(u"m", rise)) : Float64(rise)
    rise_m > 0 || throw(DomainError(rise, "VaultShell rise must be positive."))
    
    # Normalize span axis
    ax_len = sqrt(span_axis[1]^2 + span_axis[2]^2)
    ax_len > 0 || throw(DomainError(span_axis, "VaultShell span axis must be nonzero."))
    ax_dir = (span_axis[1] / ax_len, span_axis[2] / ax_len)
    perp_dir = (-ax_dir[2], ax_dir[1])
    
    # Extract corner positions (2D for triangulation)
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    z0 = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
    
    # Centroid
    cx = sum(p[1] for p in boundary_pts) / nc
    cy = sum(p[2] for p in boundary_pts) / nc
    
    # Project corners onto span axis to get span bounds
    span_projs = [ax_dir[1] * (p[1] - cx) + ax_dir[2] * (p[2] - cy) for p in boundary_pts]
    span_min, span_max = extrema(span_projs)
    L = span_max - span_min
    L > 0 || throw(DomainError(
        span_axis, "VaultShell corners have zero extent along the span axis."))
    
    # Parabolic profile: z(x,y) = z₀ + (4h/L²) · s · (L - s)
    # where s is projection of (x,y) onto span axis from span_min
    function _project_to_vault(x::Float64, y::Float64)
        s = ax_dir[1] * (x - cx) + ax_dir[2] * (y - cy) - span_min
        z = z0 + (4.0 * rise_m / L^2) * s * (L - s)
        return (x, y, z)
    end
    
    # ===========================================================================
    # Convert parameters to meters (same as Shell)
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
        5.0 * h_near
    else
        nothing
    end
    
    # Bounding box
    xmin = minimum(p[1] for p in boundary_pts)
    xmax = maximum(p[1] for p in boundary_pts)
    ymin = minimum(p[2] for p in boundary_pts)
    ymax = maximum(p[2] for p in boundary_pts)
    Lx = xmax - xmin
    Ly = ymax - ymin

    if h_far !== nothing
        _warn_mesh_density(h_far, Lx, Ly; reference_shortest_side=mesh_density_warn_shortest_side_m)
    end
    
    # ===========================================================================
    # Resolve refinement targets
    # ===========================================================================
    interior_node_pts = Tuple{Float64, Float64}[(ustrip(u"m", node.position[1]), 
                                                  ustrip(u"m", node.position[2])) 
                                                 for node in interior_nodes]
    
    refine_pts = if h_near !== nothing
        if refinement_targets isa Vector
            Tuple{Float64, Float64}[(ustrip(u"m", nd.position[1]), 
                                     ustrip(u"m", nd.position[2])) for nd in refinement_targets]
        elseif refinement_targets === :none
            Tuple{Float64, Float64}[]
        else
            # :auto — use interior_nodes as refinement targets
            interior_node_pts
        end
    else
        Tuple{Float64, Float64}[]
    end
    
    # Effective n
    effective_n = if h_far !== nothing
        max(1, ceil(Int, max(Lx, Ly) / h_far))
    else
        n
    end
    
    # ===========================================================================
    # Delaunay triangulation with optional Ruppert refinement
    # ===========================================================================
    want_refinement = h_near !== nothing && !isempty(refine_pts)
    
    tol_dedup = 1e-6
    round_coord(x) = round(Int64, x / tol_dedup)
    
    # Build node map to track corner and interior nodes
    corner_map = Dict{Tuple{Int64, Int64}, Node}()
    for corner in corners
        key = (round_coord(ustrip(u"m", corner.position[1])),
               round_coord(ustrip(u"m", corner.position[2])))
        corner_map[key] = corner
    end
    
    if want_refinement
        # ── Minimal seed + Ruppert refinement ──
        seed_points = copy(boundary_pts)
        existing_keys = Set{Tuple{Int64,Int64}}(
            (round_coord(p[1]), round_coord(p[2])) for p in boundary_pts)
        
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
        
        boundary_loop = vcat(collect(1:nc), [1])
        tri = DT.triangulate(seed_points; boundary_nodes=[boundary_loop], rng=_mesh_rng())
        
        far_h  = h_far !== nothing ? h_far : max(Lx, Ly) / max(effective_n, 1)
        A_far  = 0.433 * far_h^2
        A_near = 0.433 * h_near^2
        r_transition = r_trans !== nothing ? r_trans : 5.0 * h_near

        refine_grid = _PointGrid(refine_pts)
        function _needs_refine(_tri, T)
            i, j, k = DT.triangle_vertices(T)
            p, q, r = DT.get_point(_tri, i, j, k)
            A = abs(DT.triangle_area(p, q, r))
            tcx = (p[1] + q[1] + r[1]) / 3
            tcy = (p[2] + q[2] + r[2]) / 3
            d_min = _min_dist(refine_grid, tcx, tcy)
            α = clamp(d_min / r_transition, 0.0, 1.0)
            A_limit = A_near + α * (A_far - A_near)
            return A > A_limit
        end

        max_pts = max(5000, 20 * length(seed_points))
        try
            DT.refine!(tri; min_angle=30.0, custom_constraint=_needs_refine,
                       max_points=max_pts, rng=_mesh_rng())
        catch e
            @warn "VaultShell Ruppert refinement failed — falling back to grid mesh" exception=(e, catch_backtrace())
            tri = nothing
        end
        
        if tri !== nothing
            return _create_vault_shells_from_tri(tri, corners, interior_nodes, interior_node_map,
                                                  boundary_pts, corner_map, section, id,
                                                  edge_support_type, _project_to_vault)
        end
    end
    
    # ── Fallback: Dense grid mesh ──
    all_points, _, interior_node_map = _generate_mesh_points_with_supports(
        boundary_pts, effective_n, _SupportLine[], interior_node_pts)
    
    boundary_loop = vcat(collect(1:nc), [1])
    tri = DT.triangulate(all_points; boundary_nodes=[boundary_loop], rng=_mesh_rng())
    
    return _create_vault_shells_from_tri(tri, corners, interior_nodes, interior_node_map,
                                          boundary_pts, corner_map, section, id,
                                          edge_support_type, _project_to_vault)
end

"""
Create vault shell elements from a 2D triangulation by projecting vertices onto the parabolic surface.
"""
function _create_vault_shells_from_tri(
    tri,
    corners::Union{NTuple{4, Node}, Vector{Node}},
    interior_nodes::Vector{Node},
    interior_node_map::Dict{Int, Int},
    boundary_pts::Vector{Tuple{Float64, Float64}},
    corner_map::Dict{Tuple{Int64, Int64}, Node},
    section::ShellSection,
    id::Symbol,
    edge_support_type::Union{Symbol, Vector{Bool}},
    project_fn::Function
)
    nc = length(corners)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)
    
    points = DT.get_points(tri)
    node_map = Dict{Int, Node}()
    edge_node_indices = Set{Int}()
    
    for (idx, pt) in enumerate(points)
        if _is_on_boundary_edge(pt, boundary_pts, tol)
            push!(edge_node_indices, idx)
        end
    end
    
    # Map corner indices
    for (i, corner) in enumerate(corners)
        node_map[i] = corner
    end
    
    # Map interior nodes (use actual Node objects for structural connectivity)
    for (input_idx, point_idx) in interior_node_map
        node_map[point_idx] = interior_nodes[input_idx]
    end
    
    # Create nodes for remaining points, projected onto vault surface
    for i in (nc+1):length(points)
        haskey(node_map, i) && continue
        pt = points[i]
        
        # Check if this matches a corner (due to floating point)
        key = (round_coord(pt[1]), round_coord(pt[2]))
        if haskey(corner_map, key)
            node_map[i] = corner_map[key]
            continue
        end
        
        # Project to vault surface
        x, y, z = project_fn(pt[1], pt[2])
        pos = [x, y, z] .* u"m"
        
        fixity = i in edge_node_indices ? edge_support_type : :free
        node_map[i] = _create_node_with_fixity(pos, fixity, id)
    end
    
    # Create shell elements
    thickness_q = section.thickness * u"m"
    shells = ShellTri3[]
    
    for T in DT.each_solid_triangle(tri)
        i, j, k = DT.triangle_vertices(T)
        haskey(node_map, i) && haskey(node_map, j) && haskey(node_map, k) || continue
        
        push!(shells, ShellTri3(
            (node_map[i], node_map[j], node_map[k]),
            thickness_q, section.E * u"Pa", section.ν;
            ρ=section.ρ, id=id))
    end
    
    return shells
end

"""
    get_vault_mesh_data(corners, span_axis, rise; target_edge_length=0.25, n=20)

Get raw vertex and face data for a parabolic vault mesh (for visualization/serialization).

Uses 2D Delaunay triangulation projected onto the parabolic surface, consistent
with `VaultShell()`.

Returns a NamedTuple with:
- `vertices::Vector{NTuple{3, Float64}}`: (x, y, z) coordinates in meters
- `faces::Vector{NTuple{3, Int}}`: Triangle connectivity (1-indexed)

# Example
```julia
mesh_data = get_vault_mesh_data(corners, (1.0, 0.0), 0.75u"m"; target_edge_length=0.15)
```
"""
function get_vault_mesh_data(
    corners::Union{NTuple{4, Node}, Vector{Node}},
    span_axis::Tuple{Float64, Float64},
    rise::Union{Quantity, Real};
    target_edge_length::Union{Quantity, Real, Nothing} = 0.25u"m",
    n::Int = 20
)
    nc = length(corners)
    @assert nc == 4 "get_vault_mesh_data requires exactly 4 corner nodes"
    
    rise_m = rise isa Quantity ? Float64(ustrip(u"m", rise)) : Float64(rise)
    @assert rise_m > 0 "Rise must be positive"
    
    # Normalize span axis
    ax_len = sqrt(span_axis[1]^2 + span_axis[2]^2)
    ax_dir = (span_axis[1] / ax_len, span_axis[2] / ax_len)
    
    # Extract corner positions (2D for triangulation)
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    z0 = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
    
    cx = sum(p[1] for p in boundary_pts) / nc
    cy = sum(p[2] for p in boundary_pts) / nc
    
    # Project corners onto span axis to get span bounds
    span_projs = [ax_dir[1] * (p[1] - cx) + ax_dir[2] * (p[2] - cy) for p in boundary_pts]
    span_min, span_max = extrema(span_projs)
    L = span_max - span_min
    
    # Parabolic projection function
    function _project_to_vault(x::Float64, y::Float64)
        s = ax_dir[1] * (x - cx) + ax_dir[2] * (y - cy) - span_min
        z = z0 + (4.0 * rise_m / L^2) * s * (L - s)
        return (x, y, z)
    end
    
    # Bounding box
    xmin = minimum(p[1] for p in boundary_pts)
    xmax = maximum(p[1] for p in boundary_pts)
    ymin = minimum(p[2] for p in boundary_pts)
    ymax = maximum(p[2] for p in boundary_pts)
    Lx = xmax - xmin
    Ly = ymax - ymin
    
    # Effective n
    h_far = if target_edge_length !== nothing
        target_edge_length isa Quantity ? Float64(ustrip(u"m", target_edge_length)) : Float64(target_edge_length)
    else
        nothing
    end
    
    effective_n = if h_far !== nothing
        max(1, ceil(Int, max(Lx, Ly) / h_far))
    else
        n
    end
    
    # Generate 2D mesh points
    all_points, _, _ = _generate_mesh_points_with_supports(boundary_pts, effective_n, _SupportLine[], Tuple{Float64, Float64}[])
    
    # Triangulate
    boundary_loop = vcat(collect(1:nc), [1])
    tri = DT.triangulate(all_points; boundary_nodes=[boundary_loop], rng=_mesh_rng())
    
    # Extract vertices (projected to vault surface) and faces
    points = DT.get_points(tri)
    vertices = NTuple{3, Float64}[]
    point_to_vertex_idx = Dict{Int, Int}()
    
    for (i, pt) in enumerate(points)
        x, y, z = _project_to_vault(pt[1], pt[2])
        push!(vertices, (x, y, z))
        point_to_vertex_idx[i] = length(vertices)
    end
    
    faces = NTuple{3, Int}[]
    for T in DT.each_solid_triangle(tri)
        i, j, k = DT.triangle_vertices(T)
        if haskey(point_to_vertex_idx, i) && haskey(point_to_vertex_idx, j) && haskey(point_to_vertex_idx, k)
            push!(faces, (point_to_vertex_idx[i], point_to_vertex_idx[j], point_to_vertex_idx[k]))
        end
    end
    
    return (vertices=vertices, faces=faces)
end
