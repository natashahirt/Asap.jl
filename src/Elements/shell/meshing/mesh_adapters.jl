
"""
    Shell(triangulation, section; id=:shell, node_map=nothing, z=0.0)

Create shell elements from an external triangulation (legacy API - all nodes free).
"""
function Shell(
    tri::DT.Triangulation,
    section::ShellSection;
    id::Symbol = :shell,
    node_map::Union{Nothing, Dict{Int, Node}} = nothing,
    z::Real = 0.0
)
    points = DT.get_points(tri)
    
    # Create node map if not provided
    if isnothing(node_map)
        node_map = Dict{Int, Node}()
        for (i, pt) in enumerate(points)
            pos = [pt[1], pt[2], z] .* u"m"
            node_map[i] = Node(pos, :free, id)
        end
    end
    
    return _create_shells_from_tri(tri, node_map, section, id)
end

# ============================================================================
# Shell from DT.Triangulation with support specification
# ============================================================================

"""
    Shell(tri, section, supports; z=0.0, tol=0.01, support_type=:pinned, id=:shell)

Create shell elements from a DelaunayTriangulation with specified supports.

# Arguments
- `tri::DT.Triangulation`: Triangulation from `mesh()` or DelaunayTriangulation.jl
- `section::ShellSection`: Thickness and material properties
- `supports`: Which vertices to fix:
  - `Vector{Int}`: Specific vertex indices (1-based)
  - `:xy_plane`: Vertices where z ≈ 0 (default for flat meshes)
  - `:none`: All nodes free

# Keywords
- `z::Real=0.0`: Z-coordinate for all nodes (2D triangulations are in XY plane)
- `tol::Real=0.01`: Tolerance for plane detection (meters)
- `support_type::Symbol=:pinned`: Fixity type for supported nodes
- `id::Symbol=:shell`: Element identifier

# Example
```julia
# Create mesh and pin boundary vertices
tri = mesh((n1, n2, n3, n4), 6)
shells = Shell(tri, section, [1, 2, 3, 4])  # pin corners

# Or auto-detect (all nodes at z≈0 for flat mesh)
shells = Shell(tri, section, :xy_plane)
```
"""
function Shell(
    tri::DT.Triangulation,
    section::ShellSection,
    supports::Union{Symbol, Vector{Int}};
    z::Real = 0.0,
    tol::Real = 0.01,
    support_type::Symbol = :pinned,
    id::Symbol = :shell
)
    points = DT.get_points(tri)
    n_pts = length(points)
    
    # Resolve support indices from Symbol or Vector
    support_set = _resolve_supports_2d(supports, n_pts, z, tol)
    
    # Create nodes with appropriate fixity
    node_map = Dict{Int, Node}()
    for (i, pt) in enumerate(points)
        pos = [pt[1], pt[2], z] .* u"m"
        fixity = i ∈ support_set ? support_type : :free
        node_map[i] = Node(pos, fixity, id)
    end
    
    return _create_shells_from_tri(tri, node_map, section, id)
end

"""Resolve support specification to a Set of indices for 2D triangulations."""
function _resolve_supports_2d(supports::Vector{Int}, n_pts::Int, z::Real, tol::Real)
    return Set(supports)
end

function _resolve_supports_2d(supports::Symbol, n_pts::Int, z::Real, tol::Real)
    if supports == :xy_plane
        # For 2D triangulation at z=z, pin all if z≈0
        abs(z) < tol ? Set(1:n_pts) : Set{Int}()
    elseif supports == :none
        Set{Int}()
    else
        error("Unknown support type: $supports. Use :xy_plane, :none, or Vector{Int}")
    end
end

# ============================================================================
# Shell from Meshes.jl SimpleMesh
# ============================================================================

"""
    Shell(mesh, section, supports=:xy_plane; tol=0.01, support_type=:pinned, id=:shell)

Create shell elements from a Meshes.jl SimpleMesh.

# Arguments
- `mesh::SimpleMesh`: Triangle mesh from Meshes.jl
- `section::ShellSection`: Thickness and material properties
- `supports`: Which vertices to fix (default: `:xy_plane`):
  - `Vector{Int}`: Specific vertex indices (1-based)
  - `:xy_plane`: Vertices where z ≈ 0
  - `:xz_plane`: Vertices where y ≈ 0
  - `:yz_plane`: Vertices where x ≈ 0
  - `:none`: All nodes free

# Keywords
- `tol::Real=0.01`: Tolerance for plane detection (in mesh coordinate units)
- `support_type::Symbol=:pinned`: Fixity type for supported nodes
- `id::Symbol=:shell`: Element identifier

# Example
```julia
using Meshes

# Create or load a mesh
mesh = SimpleMesh(points, connectivity)

# Auto-pin nodes at z≈0 (common for shells sitting on XY plane)
shells = Shell(mesh, section)

# Or specify explicit support indices
shells = Shell(mesh, section, boundary_indices)

# Or use another plane
shells = Shell(mesh, section, :xz_plane)  # y≈0

model = ShellModel(get_nodes(shells), shells, loads)
solve!(model)
```
"""
function Shell(
    mesh::SimpleMesh,
    section::ShellSection,
    supports::Union{Symbol, Vector{Int}} = :xy_plane;
    tol::Real = 0.01,
    support_type::Symbol = :pinned,
    id::Symbol = :shell
)
    verts = Meshes.vertices(mesh)
    n_verts = length(verts)
    
    # Resolve support indices
    support_set = _resolve_supports_3d(supports, verts, tol)
    
    # Create nodes from mesh vertices
    nodes = Node[]
    for (i, v) in enumerate(verts)
        c = coords(v)
        # Meshes.jl coords - extract x, y and handle 2D vs 3D for z
        x = _coord_to_meters(c.x)
        y = _coord_to_meters(c.y)
        z = _get_z_coord(c)
        pos = [x, y, z] .* u"m"
        fixity = i ∈ support_set ? support_type : :free
        push!(nodes, Node(pos, fixity, id))
    end
    
    # Create shells from mesh triangles
    # Use topology() to get Connectivity objects which have indices()
    shells = ShellTri3[]
    thickness_q = section.thickness * u"m"
    
    for conn in Meshes.topology(mesh)
        idx = Meshes.indices(conn)
        # Handle different index tuple formats
        i, j, k = idx[1], idx[2], idx[3]
        push!(shells, ShellTri3(
            (nodes[i], nodes[j], nodes[k]),
            thickness_q,
            section.E * u"Pa",
            section.ν;
            ρ = section.ρ,
            id = id
        ))
    end
    
    return shells
end

# ============================================================================
# Shell from raw points + triangles (most flexible interface)
# ============================================================================

"""
    Shell(points, triangles, section, supports=:xy_plane; tol=0.01, support_type=:pinned, id=:shell)

Create shell elements from raw point coordinates and triangle connectivity.

This is the most flexible interface - no mesh library required. Points can be any
collection of 3D coordinates (tuples, vectors, or matrix rows). Triangle connectivity
specifies which points form each triangle (1-based indices).

# Arguments
- `points`: Vertex coordinates in one of these formats:
  - `Vector{NTuple{3,T}}` where T is Real or Quantity
  - `Vector{Vector{T}}` with 3-element inner vectors
  - `Matrix{T}` with size (n_points, 3)
  - Coordinates are assumed in meters if unitless
- `triangles`: Triangle connectivity in one of these formats:
  - `Vector{NTuple{3,Int}}` - e.g., `[(1,2,3), (2,4,3)]`
  - `Vector{Vector{Int}}` with 3-element inner vectors
- `section::ShellSection`: Thickness and material properties
- `supports`: Which vertices to fix (default: `:xy_plane`):
  - `Vector{Int}`: Specific vertex indices (1-based)
  - `:xy_plane`: Vertices where z ≈ 0
  - `:xz_plane`: Vertices where y ≈ 0
  - `:yz_plane`: Vertices where x ≈ 0
  - `:none`: All nodes free

# Keywords
- `tol::Real=0.01`: Tolerance for plane detection (meters)
- `support_type::Symbol=:pinned`: Fixity type for supported nodes
- `id::Symbol=:shell`: Element identifier

# Examples
```julia
# From tuples (most common)
points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.2), (1.5, 1.0, 0.2)]
triangles = [(1, 2, 3), (2, 4, 3)]
shells = Shell(points, triangles, section, [1, 2])  # pin nodes 1 and 2

# From matrix (e.g., loaded from file)
pts = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.5 1.0 0.2; 1.5 1.0 0.2]  # 4×3
tris = [(1, 2, 3), (2, 4, 3)]
shells = Shell(pts, tris, section)  # auto-pin at z≈0

# With unitful coordinates
points = [(0.0u"m", 0.0u"m", 0.0u"m"), (1.0u"m", 0.0u"m", 0.0u"m"), ...]
shells = Shell(points, triangles, section, :none)

model = ShellModel(get_nodes(shells), shells, loads)
solve!(model)
```
"""
function Shell(
    points::Union{Vector{<:NTuple{3}}, Vector{<:AbstractVector}, AbstractMatrix},
    triangles::Union{Vector{<:NTuple{3,Int}}, Vector{<:AbstractVector{Int}}},
    section::ShellSection,
    supports::Union{Symbol, Vector{Int}} = :xy_plane;
    tol::Real = 0.01,
    support_type::Symbol = :pinned,
    id::Symbol = :shell
)
    # Normalize points to Vector of 3-tuples (Float64, in meters)
    pts = _normalize_points(points)
    n_pts = length(pts)
    
    # Normalize triangles to Vector of 3-tuples
    tris = _normalize_triangles(triangles)
    
    # Resolve support indices
    support_set = _resolve_supports_raw(supports, pts, tol)
    
    # Create nodes
    nodes = Node[]
    for (i, pt) in enumerate(pts)
        pos = [pt[1], pt[2], pt[3]] .* u"m"
        fixity = i ∈ support_set ? support_type : :free
        push!(nodes, Node(pos, fixity, id))
    end
    
    # Create shell elements
    shells = ShellTri3[]
    thickness_q = section.thickness * u"m"
    
    for tri in tris
        i, j, k = tri
        @assert 1 <= i <= n_pts && 1 <= j <= n_pts && 1 <= k <= n_pts "Triangle index out of bounds"
        push!(shells, ShellTri3(
            (nodes[i], nodes[j], nodes[k]),
            thickness_q,
            section.E * u"Pa",
            section.ν;
            ρ = section.ρ,
            id = id
        ))
    end
    
    return shells
end

# --- Point normalization helpers ---

"""Normalize various point formats to Vector{NTuple{3,Float64}} in meters."""
function _normalize_points(pts::Vector{<:NTuple{3}})
    return [(_to_meters(p[1]), _to_meters(p[2]), _to_meters(p[3])) for p in pts]
end

function _normalize_points(pts::Vector{<:AbstractVector})
    return [(_to_meters(p[1]), _to_meters(p[2]), _to_meters(p[3])) for p in pts]
end

function _normalize_points(pts::AbstractMatrix)
    size(pts, 2) == 3 || error("Point matrix must have 3 columns (x, y, z)")
    return [(_to_meters(pts[i,1]), _to_meters(pts[i,2]), _to_meters(pts[i,3])) for i in 1:size(pts,1)]
end

"""Convert coordinate to Float64 meters."""
_to_meters(x::Real) = Float64(x)
_to_meters(x::Quantity) = Float64(ustrip(u"m", x))

# --- Triangle normalization helpers ---

"""Normalize triangle connectivity to Vector{NTuple{3,Int}}."""
function _normalize_triangles(tris::Vector{<:NTuple{3,Int}})
    return collect(tris)
end

function _normalize_triangles(tris::Vector{<:AbstractVector{Int}})
    return [(t[1], t[2], t[3]) for t in tris]
end

# --- Support resolution for raw points ---

"""Resolve support specification to a Set of indices for raw point arrays."""
function _resolve_supports_raw(supports::Vector{Int}, pts::Vector{NTuple{3,Float64}}, tol::Real)
    return Set(supports)
end

function _resolve_supports_raw(supports::Symbol, pts::Vector{NTuple{3,Float64}}, tol::Real)
    if supports == :xy_plane
        Set(i for (i, p) in enumerate(pts) if abs(p[3]) < tol)
    elseif supports == :xz_plane
        Set(i for (i, p) in enumerate(pts) if abs(p[2]) < tol)
    elseif supports == :yz_plane
        Set(i for (i, p) in enumerate(pts) if abs(p[1]) < tol)
    elseif supports == :none
        Set{Int}()
    else
        error("Unknown support type: $supports. Use :xy_plane, :xz_plane, :yz_plane, :none, or Vector{Int}")
    end
end

"""Resolve support specification to a Set of indices for 3D meshes."""
function _resolve_supports_3d(supports::Vector{Int}, verts, tol::Real)
    return Set(supports)
end

function _resolve_supports_3d(supports::Symbol, verts, tol::Real)
    if supports == :xy_plane
        Set(i for (i, v) in enumerate(verts) if abs(_get_z(v)) < tol)
    elseif supports == :xz_plane
        Set(i for (i, v) in enumerate(verts) if abs(_get_y(v)) < tol)
    elseif supports == :yz_plane
        Set(i for (i, v) in enumerate(verts) if abs(_get_x(v)) < tol)
    elseif supports == :none
        Set{Int}()
    else
        error("Unknown support type: $supports. Use :xy_plane, :xz_plane, :yz_plane, :none, or Vector{Int}")
    end
end

# Coordinate extractors for Meshes.jl Points
_get_x(v) = _coord_to_meters(coords(v).x)
_get_y(v) = _coord_to_meters(coords(v).y)
_get_z(v) = _get_z_coord(coords(v))

"""Convert Meshes.jl coordinate (may be unitful) to Float64 meters."""
_coord_to_meters(x::Real) = Float64(x)
_coord_to_meters(x::Quantity) = Float64(ustrip(u"m", x))

"""Extract z coordinate, returning 0.0 for 2D coords."""
function _get_z_coord(c)
    # Check if this is a 3D coordinate system by trying to access z
    # Meshes.jl Cartesian2D throws BoundsError on .z access
    if hasproperty(c, :z)
        try
            return _coord_to_meters(c.z)
        catch e
            e isa BoundsError && return 0.0
            rethrow(e)
        end
    end
    return 0.0
end

# ============================================================================
# get_nodes - Extract unique nodes from shells
# ============================================================================

"""
    get_nodes(shells)

Extract all unique nodes from shell elements.

```julia
shells = Shell(corners, section)
model = ShellModel(get_nodes(shells), shells, loads)
```
"""
function get_nodes(shells::Vector{ShellTri3})
    seen = Set{UInt64}()
    nodes = Node[]
    
    for shell in shells
        for node in shell.nodes
            id = objectid(node)
            if !(id in seen)
                push!(seen, id)
                push!(nodes, node)
            end
        end
    end
    
    return nodes
end

function get_nodes(elements::Vector{<:AbstractElement})
    seen = Set{UInt64}()
    nodes = Node[]
    
    for elem in elements
        if hasproperty(elem, :nodes)
            for node in elem.nodes
                id = objectid(node)
                if !(id in seen)
                    push!(seen, id)
                    push!(nodes, node)
                end
            end
        elseif hasproperty(elem, :node_i) && hasproperty(elem, :node_j)
            for node in (elem.node_i, elem.node_j)
                id = objectid(node)
                if !(id in seen)
                    push!(seen, id)
                    push!(nodes, node)
                end
            end
        end
    end
    
    return nodes
end
