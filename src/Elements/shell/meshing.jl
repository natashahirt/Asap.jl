#=
Shell Creation Utilities
========================

Convenient functions to create shell elements from corner nodes or external meshes.
Uses DelaunayTriangulation.jl for automatic polygon triangulation.
=#

import DelaunayTriangulation as DT

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

# Common presets
const Concrete_ShellSection_150mm = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete_150mm)
const Concrete_ShellSection_200mm = ShellSection(0.20u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete_200mm)

# ============================================================================
# DiaphragmSection - rigid in-plane, massless
# ============================================================================

"""
    DiaphragmSection(; thickness=0.01u"m", E=1e12u"Pa", ν=0.3)

A section for rigid diaphragm elements: very stiff in-plane, zero mass.

Diaphragms constrain all connected nodes to move together in-plane (rigid body motion)
while not contributing mass or significant out-of-plane stiffness.

# Use Cases
- Floor diaphragms in lateral analysis
- Connecting frames at floor levels
- Distributing lateral loads to frames

# Properties
- `E`: Very high (1 TPa default) → rigid in-plane
- `ρ`: Zero → no mass contribution
- `thickness`: Small (10mm default) → minimal out-of-plane stiffness

# Example
```julia
diaphragm_sec = DiaphragmSection()
diaphragms = Diaphragm((n1, n2, n3, n4), diaphragm_sec)
```
"""
struct DiaphragmSection
    thickness::Float64  # meters (small)
    E::Float64          # Pa (very high)
    ν::Float64
    name::Symbol
    
    function DiaphragmSection(;
        thickness::Quantity = 0.01u"m",
        E::Quantity = 1e12u"Pa",  # 1 TPa - effectively rigid
        ν::Real = 0.3,
        name::Symbol = :diaphragm
    )
        @assert 0.0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        t = ustrip(u"m", thickness)
        E_val = ustrip(u"Pa", E)
        return new(t, E_val, Float64(ν), name)
    end
end

# Default diaphragm section
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
    @assert nc >= 3 "Shell requires at least 3 corner nodes"
    
    # Convert interior supports to _SupportLine
    support_lines = [_to_support_line(s) for s in interior_supports]
    
    # Extract corner positions
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    
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

Create shell elements from an external triangulation.
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
