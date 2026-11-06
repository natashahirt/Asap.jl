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
    
    all_points, _, _ = _generate_mesh_points_with_supports(boundary_pts, n, _SupportLine[])
    
    boundary_loop = vcat(collect(1:nc), [1])
    tri = DT.triangulate(all_points; boundary_nodes=[boundary_loop], rng=_mesh_rng())
    
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
# SupportLine - helper constructors and converters
# ============================================================================
# The struct itself is forward-declared above _collect_refinement_targets.

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

