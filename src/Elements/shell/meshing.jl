#=
Shell Creation Utilities
========================

Convenient functions to create shell elements from corner nodes or external meshes.
Uses DelaunayTriangulation.jl for automatic polygon triangulation.
=#

import DelaunayTriangulation as DT
import Meshes: SimpleMesh, Point, coords

# ============================================================================
# ShellSection - combines thickness and material properties
# ============================================================================

"""
    ShellSection(thickness, E, ν; ρ=0.0u"kg/m^3", name=:section)
    ShellSection(thickness, material::ShellMaterial)

A shell section combining thickness with material properties.

# Examples
```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
section = ShellSection(0.15u"m", concrete_material)
```
"""
struct ShellSection
    thickness::Float64  # meters
    E::Float64          # Pa
    ν::Float64
    ρ::Float64          # kg/m³
    name::Symbol
    
    function ShellSection(
        thickness::Quantity,
        E::Quantity,
        ν::Real;
        ρ::Quantity = 0.0u"kg/m^3",
        name::Symbol = :section
    )
        @assert 0.0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        t = ustrip(u"m", thickness)
        E_val = ustrip(u"Pa", E)
        ρ_val = ustrip(u"kg/m^3", ρ)
        return new(t, E_val, Float64(ν), ρ_val, name)
    end
    
    function ShellSection(thickness::Quantity, material::ShellMaterial)
        t = ustrip(u"m", thickness)
        E_val = ustrip(material.E)
        ρ_val = ustrip(material.ρ)
        return new(t, E_val, material.ν, ρ_val, material.name)
    end
end

# Note: ShellSection presets removed. Create sections inline:
# ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete)

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

# ============================================================================
# DiaphragmSection - rigid in-plane, massless
# ============================================================================

"""
    DiaphragmSection(; thickness=0.01u"m", E=1e12u"Pa", ν=0.3, name=:diaphragm)

A section for rigid diaphragm elements: very stiff in-plane, zero mass.

Diaphragms constrain all connected nodes to move together in-plane (rigid body motion)
while not contributing mass or significant out-of-plane stiffness.

# Use Cases
- Floor diaphragms in lateral analysis
- Connecting frames at floor levels
- Distributing lateral loads to frames

# Design Parameters

**thickness** (default: 10mm = 0.01m): Controls out-of-plane stiffness.
- Smaller thickness → more flexible out-of-plane (desirable for pure in-plane constraint)
- Typical range: 5-20mm
- Note: Too thin can cause numerical issues; too thick adds unwanted bending stiffness

**E** (default: 1 TPa = 1e12 Pa): Young's modulus for in-plane stiffness.
- Higher E → more rigid in-plane behavior
- 1 TPa is ~5× steel stiffness, sufficient for practical rigidity
- Use 200 GPa for steel-like diaphragm, 30 GPa for concrete-like

**ν** (default: 0.3): Poisson's ratio
- Affects in-plane lateral contraction
- 0.3 is typical for metals; 0.2 for concrete

# Example
```julia
# Default rigid diaphragm
diaphragm_sec = DiaphragmSection()

# Steel deck diaphragm (less rigid, for sensitivity studies)
steel_diaphragm = DiaphragmSection(E=200u"GPa", thickness=0.003u"m")

# Create diaphragm elements
diaphragms = Diaphragm((n1, n2, n3, n4), diaphragm_sec)
```

# Notes
- Mass (ρ) is always 0 for diaphragms - they don't contribute to modal analysis
- If you need mass, use `Shell()` with full material properties instead
"""
struct DiaphragmSection
    thickness::Float64  # meters
    E::Float64          # Pa
    ν::Float64
    name::Symbol
    
    function DiaphragmSection(;
        thickness::Quantity = 0.01u"m",
        E::Quantity = 1e12u"Pa",
        ν::Real = 0.3,
        name::Symbol = :diaphragm
    )
        @assert 0.0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        t = ustrip(u"m", thickness)
        E_val = ustrip(u"Pa", E)
        return new(t, E_val, Float64(ν), name)
    end
end

# Default diaphragm section (1 TPa, 10mm thick)
const RigidDiaphragm = DiaphragmSection()

# ============================================================================
# Diaphragm() - Create rigid diaphragm elements
# ============================================================================

"""
    Diaphragm(corners, section; n=2, id=:diaphragm)
    Diaphragm(corners; n=2, id=:diaphragm)

Create rigid diaphragm elements from corner nodes.

Diaphragms are shell elements with:
- **Very high in-plane stiffness** → all nodes move together in XY
- **Zero mass** → no contribution to dynamic analysis mass matrix
- **Minimal out-of-plane stiffness** → doesn't affect vertical behavior

# Arguments
- `corners`: Tuple or vector of corner `Node`s (3+ nodes)
- `section`: `DiaphragmSection` (optional, defaults to `RigidDiaphragm`)
- `n`: Mesh refinement (default: 2, coarser than Shell since we just need connectivity)

# Example
```julia
# Floor diaphragm connecting frames
floor_corners = (col1_top, col2_top, col3_top, col4_top)
diaphragms = Diaphragm(floor_corners)

# With custom stiffness
soft_diaphragm = DiaphragmSection(E=200u"GPa")  # Steel-like
diaphragms = Diaphragm(corners, soft_diaphragm)
```

# When to Use
- **Lateral analysis**: Connect frames at floor levels
- **Load distribution**: Ensure lateral loads go to all frames
- **Simplified models**: When slab moments aren't needed

# When NOT to Use
- When you need slab bending moments → use `Shell()`
- When mass matters for dynamics → use `Shell()` with proper ρ
"""
function Diaphragm(
    corners::Union{NTuple{N, Node}, Vector{Node}},
    section::DiaphragmSection;
    n::Int = 2,
    id::Symbol = :diaphragm
) where N
    nc = length(corners)
    @assert nc >= 3 "Diaphragm requires at least 3 corner nodes"
    
    # Use the Shell meshing infrastructure, but with diaphragm properties
    # Get triangulation
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    
    all_points, _ = _generate_mesh_points_with_supports(boundary_pts, n, _SupportLine[])
    
    boundary_loop = vcat(collect(1:nc), [1])
    tri = DT.triangulate(all_points; boundary_nodes=[boundary_loop])
    
    # Compute average z
    z_val = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
    
    # Create node map
    points = DT.get_points(tri)
    node_map = Dict{Int, Node}()
    
    # Map corner nodes
    for (i, corner) in enumerate(corners)
        node_map[i] = corner
    end
    
    # Create interior nodes (free, since diaphragm doesn't need edge supports)
    for i in (nc+1):length(points)
        pt = points[i]
        pos = [pt[1], pt[2], z_val] .* u"m"
        node_map[i] = Node(pos, :free, id)
    end
    
    # Create diaphragm shell elements with zero density
    diaphragms = ShellTri3[]
    thickness_q = section.thickness * u"m"
    
    for triangle in DT.each_solid_triangle(tri)
        i, j, k = triangle
        push!(diaphragms, ShellTri3(
            (node_map[i], node_map[j], node_map[k]), 
            thickness_q, 
            section.E * u"Pa", 
            section.ν; 
            ρ = 0.0,  # Zero mass!
            id = id
        ))
    end
    
    return diaphragms
end

# Default section version
function Diaphragm(
    corners::Union{NTuple{N, Node}, Vector{Node}};
    n::Int = 2,
    id::Symbol = :diaphragm
) where N
    return Diaphragm(corners, RigidDiaphragm; n=n, id=id)
end

# ============================================================================
# SupportLine - internal representation of a support line
# ============================================================================

"""
Internal struct representing a support line for meshing.
"""
struct _SupportLine
    start_pt::Tuple{Float64, Float64}  # (x, y) in meters
    end_pt::Tuple{Float64, Float64}
end

"""
Extract support line geometry from an Element.
"""
function _support_line_from_element(elem::Element)
    p1 = elem.nodeStart.position
    p2 = elem.nodeEnd.position
    start_pt = (ustrip(u"m", p1[1]), ustrip(u"m", p1[2]))
    end_pt = (ustrip(u"m", p2[1]), ustrip(u"m", p2[2]))
    return _SupportLine(start_pt, end_pt)
end

"""
Extract support line geometry from a node pair.
"""
function _support_line_from_nodes(n1::Node, n2::Node)
    p1, p2 = n1.position, n2.position
    start_pt = (ustrip(u"m", p1[1]), ustrip(u"m", p1[2]))
    end_pt = (ustrip(u"m", p2[1]), ustrip(u"m", p2[2]))
    return _SupportLine(start_pt, end_pt)
end

"""
Convert various input types to _SupportLine.
"""
function _to_support_line(input)
    if input isa Element
        return _support_line_from_element(input)
    elseif input isa Tuple{Node, Node}
        return _support_line_from_nodes(input[1], input[2])
    elseif input isa _SupportLine
        return input
    else
        error("interior_supports must contain Elements or (Node, Node) tuples")
    end
end

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
    
    # Convert interior supports to _SupportLine
    support_lines = [_to_support_line(s) for s in interior_supports]
    
    # Generate all points (boundary + interior + support lines) without duplicates
    all_points, support_point_indices = _generate_mesh_points_with_supports(boundary_pts, n, support_lines)
    
    # Build constrained triangulation
    boundary_loop = vcat(collect(1:nc), [1])
    return DT.triangulate(all_points; boundary_nodes=[boundary_loop])
end

# ============================================================================
# Shell() - Create shell elements
# ============================================================================

"""
    Shell(corners, section; n=4, id=:shell, interior_supports=[], edge_support_type=:pinned, interior_support_type=:pinned)
    Shell(corners, n, section; ...)

Create triangular shell elements from any polygon with support conditions.

# Arguments
- `corners`: Tuple or vector of corner `Node`s in counter-clockwise order
- `section`: `ShellSection` defining thickness and material properties
- `n`: Subdivision level (default: 4). Higher = finer mesh.

# Keyword Arguments
- `id::Symbol`: Element identifier (default: `:shell`)
- `interior_supports`: List of interior support lines - Elements or (Node, Node) pairs
- `edge_support_type::Symbol`: Fixity for edge nodes (default: `:pinned`)
- `interior_support_type::Symbol`: Fixity for interior support line nodes (default: `:pinned`)

# Support Types
- `:pinned` - x, y, z translations fixed, rotations free (simply-supported)
- `:zfixed` - only z translation fixed (vertical support, free to slide)
- `:fixed` - all 6 DOFs fixed (clamped)
- `:free` - no constraints (for edges that aren't supported)
- `Vector{Bool}` - custom DOF pattern [x, y, z, θx, θy, θz]

# Example
```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2)

# Simple slab (edges supported)
shells = Shell((n1, n2, n3, n4), section)

# Multi-bay slab with interior beam
shells = Shell((n1, n2, n3, n4), section;
    interior_supports = [interior_beam],
    edge_support_type = :pinned,
    interior_support_type = :pinned
)

model = ShellModel(get_nodes(shells), shells, loads)
```
"""
function Shell(
    corners::Union{NTuple{N, Node}, Vector{Node}},
    n::Int,
    section::ShellSection;
    id::Symbol = :shell,
    interior_supports::Vector = [],
    edge_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    interior_support_type::Union{Symbol, Vector{Bool}} = :pinned
) where N
    nc = length(corners)
    
    # ===========================================================================
    # Input Validation
    # ===========================================================================
    @assert nc >= 3 "Shell requires at least 3 corner nodes, got $nc"
    @assert n >= 1 "Mesh refinement n must be at least 1, got $n"
    
    # Extract corner positions for validation
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    
    # Check for degenerate polygon (zero or negative area)
    area = _shoelace_area(boundary_pts)
    @assert abs(area) > 1e-12 "Shell polygon is degenerate (zero area). Check corner positions."
    
    # Check for self-intersection (simple polygon test)
    @assert !_is_self_intersecting(boundary_pts) "Shell polygon is self-intersecting. Corners must form a simple polygon."
    
    # ===========================================================================
    # Mesh Generation
    # ===========================================================================
    
    # Convert interior supports to _SupportLine
    support_lines = [_to_support_line(s) for s in interior_supports]
    
    # Generate mesh points with supports
    all_points, support_point_indices = _generate_mesh_points_with_supports(boundary_pts, n, support_lines)
    
    # Build constrained triangulation
    boundary_loop = vcat(collect(1:nc), [1])
    tri = DT.triangulate(all_points; boundary_nodes=[boundary_loop])
    
    # Compute average z for interior nodes
    z_val = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
    
    # Create node map
    points = DT.get_points(tri)
    node_map = Dict{Int, Node}()
    edge_node_indices = Set{Int}()
    
    # Track which point indices are on edges
    tol = 1e-6
    for (idx, pt) in enumerate(points)
        if _is_on_boundary_edge(pt, boundary_pts, tol)
            push!(edge_node_indices, idx)
        end
    end
    
    # Map corner nodes (indices 1:nc) - corners keep their original fixity
    for (i, corner) in enumerate(corners)
        node_map[i] = corner
    end
    
    # Create non-corner nodes
    for i in (nc+1):length(points)
        pt = points[i]
        pos = [pt[1], pt[2], z_val] .* u"m"
        
        # Determine fixity based on location
        if i in support_point_indices
            # Interior support line node
            fixity = interior_support_type
        elseif i in edge_node_indices
            # Edge node (not corner)
            fixity = edge_support_type
        else
            # Pure interior node
            fixity = :free
        end
        
        node_map[i] = _create_node_with_fixity(pos, fixity, id)
    end
    
    # Create shell elements
    return _create_shells_from_tri(tri, node_map, section, id)
end

# Default n=4
function Shell(
    corners::Union{NTuple{N, Node}, Vector{Node}},
    section::ShellSection;
    n::Int = 4,
    id::Symbol = :shell,
    interior_supports::Vector = [],
    edge_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    interior_support_type::Union{Symbol, Vector{Bool}} = :pinned
) where N
    return Shell(corners, n, section; 
        id=id, 
        interior_supports=interior_supports,
        edge_support_type=edge_support_type,
        interior_support_type=interior_support_type
    )
end

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

# ============================================================================
# Internal helpers
# ============================================================================

"""
Create a Node with the specified fixity.
"""
function _create_node_with_fixity(pos, fixity::Symbol, id::Symbol)
    return Node(pos, fixity, id)
end

function _create_node_with_fixity(pos, fixity::Vector{Bool}, id::Symbol)
    @assert length(fixity) == 6 "Custom fixity must be a 6-element Bool vector"
    node = Node(pos, :free, id)
    node.dof .= fixity
    return node
end

"""
Check if a point is on any boundary edge.
"""
function _is_on_boundary_edge(pt::Tuple{Float64, Float64}, boundary_pts::Vector{Tuple{Float64, Float64}}, tol::Float64)
    nc = length(boundary_pts)
    for i in 1:nc
        p1 = boundary_pts[i]
        p2 = boundary_pts[mod1(i+1, nc)]
        if _point_on_segment(pt, p1, p2, tol)
            return true
        end
    end
    return false
end

"""
Check if a point lies on a line segment.
"""
function _point_on_segment(pt::Tuple{Float64, Float64}, p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, tol::Float64)
    # Check if pt is collinear with p1-p2 and within the segment bounds
    dx, dy = p2[1] - p1[1], p2[2] - p1[2]
    len = hypot(dx, dy)
    if len < tol
        return hypot(pt[1] - p1[1], pt[2] - p1[2]) < tol
    end
    
    # Project pt onto line
    t = ((pt[1] - p1[1]) * dx + (pt[2] - p1[2]) * dy) / (len * len)
    
    # Check if within segment
    if t < -tol/len || t > 1 + tol/len
        return false
    end
    
    # Check distance from line
    proj_x = p1[1] + t * dx
    proj_y = p1[2] + t * dy
    dist = hypot(pt[1] - proj_x, pt[2] - proj_y)
    
    return dist < tol
end

"""
Create shells from triangulation and node map.
"""
function _create_shells_from_tri(
    tri::DT.Triangulation,
    node_map::Dict{Int, Node},
    section::ShellSection,
    id::Symbol
)
    shells = ShellTri3[]
    thickness_q = section.thickness * u"m"
    
    for triangle in DT.each_solid_triangle(tri)
        i, j, k = triangle
        push!(shells, ShellTri3(
            (node_map[i], node_map[j], node_map[k]), 
            thickness_q, 
            section.E * u"Pa", 
            section.ν; 
            ρ=section.ρ,
            id=id
        ))
    end
    
    return shells
end

"""
Generate mesh points including support line points.
Returns (all_points, support_point_indices).
"""
function _generate_mesh_points_with_supports(
    boundary_pts::Vector{Tuple{Float64, Float64}}, 
    n::Int,
    support_lines::Vector{_SupportLine}
)
    nc = length(boundary_pts)
    
    # Tolerance for deduplication (in meters)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)
    
    # Use a dict to dedupe: key is rounded coords, value is (exact coords, is_support)
    points_dict = Dict{Tuple{Int64, Int64}, Tuple{Tuple{Float64, Float64}, Bool}}()
    
    # Add boundary points first (these take priority, not supports)
    for pt in boundary_pts
        key = (round_coord(pt[1]), round_coord(pt[2]))
        points_dict[key] = (pt, false)
    end
    
    # Add support line points
    for line in support_lines
        line_pts = _points_along_line(line.start_pt, line.end_pt, n)
        for pt in line_pts
            key = (round_coord(pt[1]), round_coord(pt[2]))
            if !haskey(points_dict, key)
                # Only add if inside polygon
                if _point_inside_polygon(pt, boundary_pts) || _is_on_boundary_edge(pt, boundary_pts, tol)
                    points_dict[key] = (pt, true)  # Mark as support point
                end
            end
        end
    end
    
    if n <= 1
        # Build result with support tracking
        result = collect(boundary_pts)
        support_indices = Set{Int}()
        
        for (key, (pt, is_support)) in points_dict
            is_boundary = any(abs(pt[1] - bp[1]) < tol && abs(pt[2] - bp[2]) < tol for bp in boundary_pts)
            if !is_boundary
                push!(result, pt)
                if is_support
                    push!(support_indices, length(result))
                end
            end
        end
        
        return result, support_indices
    end
    
    # Compute bounding box
    xmin = minimum(p[1] for p in boundary_pts)
    xmax = maximum(p[1] for p in boundary_pts)
    ymin = minimum(p[2] for p in boundary_pts)
    ymax = maximum(p[2] for p in boundary_pts)
    
    dx = (xmax - xmin) / n
    dy = (ymax - ymin) / n
    
    # Add edge subdivision points (not supports)
    for i in 1:nc
        p1 = boundary_pts[i]
        p2 = boundary_pts[mod1(i+1, nc)]
        
        for k in 1:n-1
            t = k / n
            pt = (p1[1] * (1-t) + p2[1] * t, p1[2] * (1-t) + p2[2] * t)
            key = (round_coord(pt[1]), round_coord(pt[2]))
            if !haskey(points_dict, key)
                points_dict[key] = (pt, false)
            end
        end
    end
    
    # Add interior grid points (not supports)
    for j in 1:n-1
        for i in 1:n-1
            x = xmin + i * dx
            y = ymin + j * dy
            pt = (x, y)
            
            if _point_inside_polygon(pt, boundary_pts)
                key = (round_coord(pt[1]), round_coord(pt[2]))
                if !haskey(points_dict, key)
                    points_dict[key] = (pt, false)
                end
            end
        end
    end
    
    # Build result: boundary points first, then others
    result = collect(boundary_pts)
    support_indices = Set{Int}()
    
    for (key, (pt, is_support)) in points_dict
        is_boundary = any(abs(pt[1] - bp[1]) < tol && abs(pt[2] - bp[2]) < tol for bp in boundary_pts)
        if !is_boundary
            push!(result, pt)
            if is_support
                push!(support_indices, length(result))
            end
        end
    end
    
    return result, support_indices
end

"""
Generate points along a line at spacing consistent with n divisions.
"""
function _points_along_line(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, n::Int)
    points = Tuple{Float64, Float64}[]
    
    # Include endpoints and n-1 interior points
    for k in 0:n
        t = k / n
        pt = (p1[1] * (1-t) + p2[1] * t, p1[2] * (1-t) + p2[2] * t)
        push!(points, pt)
    end
    
    return points
end

"""
Check if point is inside polygon (ray casting).
"""
function _point_inside_polygon(pt::Tuple{Float64, Float64}, polygon::Vector{Tuple{Float64, Float64}})
    x, y = pt
    n = length(polygon)
    inside = false
    
    j = n
    for i in 1:n
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    
    return inside
end
