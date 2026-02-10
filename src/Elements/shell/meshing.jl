#=
Shell Creation Utilities
========================

Convenient functions to create shell elements from corner nodes or external meshes.
Uses DelaunayTriangulation.jl for automatic polygon triangulation.
Includes structured mesh generation (T3block) adapted from FinEtools.jl (Petr Krysl, MIT License).
=#

import DelaunayTriangulation as DT
import Meshes
import Meshes: SimpleMesh, Point, coords

# ============================================================================
# Structured Mesh Generation (adapted from FinEtools.jl MeshTriangleModule)
# ============================================================================

"""
    _t3block(Lx, Ly, nx, ny) -> (points, connectivity)

Generate a structured triangular mesh of a rectangle [0,Lx] × [0,Ly].

Returns `(points, conn)` where:
- `points::Vector{Tuple{Float64,Float64}}` — node coordinates
- `conn::Vector{NTuple{3,Int}}` — triangle connectivity (1-indexed)

Each rectangular cell is split into 2 triangles (lower-left diagonal).
Produces `(nx+1)*(ny+1)` nodes and `2*nx*ny` triangles.

Adapted from FinEtools.jl `T3block` by Petr Krysl (MIT License).
"""
function _t3block(Lx::Float64, Ly::Float64, nx::Int, ny::Int)
    xs = range(0.0, Lx, length=nx+1)
    ys = range(0.0, Ly, length=ny+1)
    return _t3blockx(collect(xs), collect(ys))
end

"""
    _t3blockx(xs, ys) -> (points, connectivity)

Generate a structured triangular mesh from explicit coordinate arrays.
Supports non-uniform (graded) spacing.

**Alternating diagonals** (checkerboard pattern): each quad cell `(i,j)` is
split along the `/` diagonal when `(i+j)` is even and along the `\\` diagonal
when `(i+j)` is odd. This eliminates the directional stiffness bias that a
uniform diagonal orientation would introduce in FE analysis.

Returns `(points, conn)` where:
- `points::Vector{Tuple{Float64,Float64}}` — node coordinates
- `conn::Vector{NTuple{3,Int}}` — triangle connectivity (1-indexed)

Adapted from FinEtools.jl `T3blockx` by Petr Krysl (MIT License), with
alternating-diagonal enhancement for structural analysis.
"""
function _t3blockx(xs::Vector{Float64}, ys::Vector{Float64})
    nL = length(xs) - 1
    nW = length(ys) - 1
    nnodes = (nL + 1) * (nW + 1)
    ncells = 2 * nL * nW

    points = Vector{Tuple{Float64, Float64}}(undef, nnodes)
    conn = Vector{NTuple{3, Int}}(undef, ncells)

    # Generate node coordinates
    f = 1
    for j in 1:(nW+1)
        for i in 1:(nL+1)
            points[f] = (xs[i], ys[j])
            f += 1
        end
    end

    # Generate triangle connectivity with alternating diagonals
    # Node numbering for cell (i, j):
    #   nw = f+(nL+1)   ne = f+(nL+1)+1
    #   sw = f           se = f+1
    gc = 1
    for i in 1:nL
        for j in 1:nW
            sw = (j - 1) * (nL + 1) + i
            se = sw + 1
            nw = sw + (nL + 1)
            ne = nw + 1
            
            if iseven(i + j)
                # "/" diagonal: sw → ne
                conn[gc]     = (sw, se, ne)
                conn[gc + 1] = (sw, ne, nw)
            else
                # "\" diagonal: se → nw
                conn[gc]     = (sw, se, nw)
                conn[gc + 1] = (se, ne, nw)
            end
            gc += 2
        end
    end

    return points, conn
end

# ============================================================================
# Graded Spacing for Mesh Refinement
# ============================================================================

"""
    _graded_spacing(x0, x1, targets, h_far, h_near, r_transition) -> Vector{Float64}

Build a non-uniform coordinate array that clusters points near `targets`.
Spacing is **symmetric** around each target — the algorithm walks from both ends
of every inter-target segment toward the midpoint and merges the results.

- `x0, x1`: domain bounds
- `targets`: locations requiring refinement (e.g., column coordinates)
- `h_far`: far-field element size
- `h_near`: near-target element size (must be ≤ h_far)
- `r_transition`: distance over which spacing grades from h_near to h_far
"""
function _graded_spacing(x0::Float64, x1::Float64, targets::Vector{Float64},
                         h_far::Float64, h_near::Float64, r_transition::Float64)
    @assert h_near <= h_far "h_near ($h_near) must be ≤ h_far ($h_far)"
    @assert x0 < x1 "x0 must be < x1"
    @assert r_transition > 0 "r_transition must be > 0"
    
    # Protected target points (must survive dedup)
    target_set = Set{Float64}([x0, x1])
    for t in targets
        if x0 ≤ t ≤ x1
            push!(target_set, t)
        end
    end
    
    # Sorted list of mandatory "waypoints" (endpoints + interior targets)
    waypoints = sort!(collect(target_set))
    
    # Collect all points
    pts = Set{Float64}()
    for wp in waypoints
        push!(pts, wp)
    end
    
    # Helper: local element size at position x
    _h_at(x) = begin
        d_min = Inf
        for t in waypoints
            d_min = min(d_min, abs(x - t))
        end
        ratio = clamp(d_min / r_transition, 0.0, 1.0)
        h_near + ratio * (h_far - h_near)
    end
    
    # For each segment [a, b], walk from BOTH ends toward the midpoint.
    # This ensures spacing is symmetric around each waypoint.
    for seg in 1:length(waypoints)-1
        a, b = waypoints[seg], waypoints[seg+1]
        mid = (a + b) / 2.0
        
        # Forward walk: a → mid
        x = a
        while x < mid - h_near * 0.05
            h_local = _h_at(x)
            x_next = min(x + h_local, mid)
            push!(pts, x_next)
            x = x_next
        end
        
        # Backward walk: b → mid
        x = b
        while x > mid + h_near * 0.05
            h_local = _h_at(x)
            x_next = max(x - h_local, mid)
            push!(pts, x_next)
            x = x_next
        end
    end
    
    result = sort!(collect(pts))
    
    # Dedup close points, but never remove protected targets
    min_spacing = h_near * 0.3
    filtered = Float64[result[1]]
    for i in 2:length(result)
        gap = result[i] - filtered[end]
        is_protected = result[i] in target_set
        if gap > min_spacing || i == length(result) || is_protected
            push!(filtered, result[i])
        end
    end
    if filtered[end] != x1
        push!(filtered, x1)
    end
    
    return filtered
end

"""
    _warn_mesh_density(h, Lx, Ly)

Warn if the effective mesh density seems unreasonably coarse or fine.
Thresholds are expressed in terms of effective `n` (divisions along the shortest side).
"""
function _warn_mesh_density(h::Float64, Lx::Float64, Ly::Float64)
    L_min = min(Lx, Ly)
    effective_n = L_min / h
    if effective_n < 4
        @warn "Shell mesh is very coarse (effective n ≈ $(round(effective_n, digits=1)) " *
              "on shortest side). Consider reducing target_edge_length."
    elseif effective_n > 100
        @warn "Shell mesh is very fine (effective n ≈ $(round(effective_n, digits=1)) " *
              "on shortest side). Consider increasing target_edge_length for faster analysis."
    end
end

"""
    _uniform_spacing(x0, x1, h) -> Vector{Float64}

Build a uniform coordinate array with target spacing `h`.
"""
function _uniform_spacing(x0::Float64, x1::Float64, h::Float64)
    n = max(1, ceil(Int, (x1 - x0) / h))
    return collect(range(x0, x1, length=n+1))
end

"""
    _is_rectangular(boundary_pts; tol=1e-6) -> Bool

Check if a polygon is a (possibly rotated) axis-aligned rectangle.
"""
function _is_rectangular(boundary_pts::Vector{Tuple{Float64, Float64}}; tol::Float64=1e-6)
    length(boundary_pts) != 4 && return false
    
    # Check that opposite edges are parallel and equal length
    p1, p2, p3, p4 = boundary_pts
    
    # Edge vectors
    e1 = (p2[1] - p1[1], p2[2] - p1[2])
    e2 = (p3[1] - p2[1], p3[2] - p2[2])
    e3 = (p4[1] - p3[1], p4[2] - p3[2])
    e4 = (p1[1] - p4[1], p1[2] - p4[2])
    
    # Check axis-aligned: each edge should be horizontal or vertical
    h1 = abs(e1[2]) < tol  # horizontal
    v1 = abs(e1[1]) < tol  # vertical
    h2 = abs(e2[2]) < tol
    v2 = abs(e2[1]) < tol
    h3 = abs(e3[2]) < tol
    v3 = abs(e3[1]) < tol
    h4 = abs(e4[2]) < tol
    v4 = abs(e4[1]) < tol
    
    return (h1 || v1) && (h2 || v2) && (h3 || v3) && (h4 || v4)
end

# ============================================================================
# SupportLine - internal representation of a support line (forward declaration)
# ============================================================================
# Defined here because _collect_refinement_targets (below) references the type.
# Helper constructors and converters follow later in the file.

"""
Internal struct representing a support line for meshing.
"""
struct _SupportLine
    start_pt::Tuple{Float64, Float64}  # (x, y) in meters
    end_pt::Tuple{Float64, Float64}
end

"""
    _collect_refinement_targets(interior_nodes, support_lines, patches)
        -> Vector{Tuple{Float64,Float64}}

Auto-detect mesh refinement targets from interior features only.
Returns point locations (in meters) where the mesh should be refined.

Targets are:
1. Interior nodes (column connection points)
2. Support-line endpoints (beam connections)
3. Patch centroids (column footprints — stress concentrations)

Boundary conditions (`edge_support_type`) are **not** considered here.
Mesh refinement and DOF fixity are orthogonal concerns.
"""
function _collect_refinement_targets(
    interior_nodes::Vector{Node},
    support_lines::Vector{_SupportLine},
    patches::Vector  # Vector{ShellPatch} — untyped to avoid forward-declaration order
)
    targets = Tuple{Float64, Float64}[]
    
    # 1. Interior nodes (column connection points)
    for node in interior_nodes
        push!(targets, (ustrip(u"m", node.position[1]), 
                        ustrip(u"m", node.position[2])))
    end
    
    # 2. Support line endpoints (beam connections)
    for line in support_lines
        push!(targets, line.start_pt)
        push!(targets, line.end_pt)
    end
    
    # 3. Patch centroids (column footprints — stress concentrations)
    for patch in patches
        push!(targets, patch.center)
    end
    
    # Deduplicate
    return unique(targets)
end

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
    
    all_points, _, _ = _generate_mesh_points_with_supports(boundary_pts, n, _SupportLine[])
    
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

# ============================================================================
# ShellPatch - stiffened region within a shell mesh
# ============================================================================

"""
    ShellPatch(polygon, section; id=:patch)
    ShellPatch(cx, cy, c1, c2, section; id=:patch)
    ShellPatch(vertices, center, section, id)

A stiffened region within a shell mesh, typically representing a column footprint.

The polygon perimeter becomes **constrained edges** in the Delaunay triangulation,
ensuring element boundaries align with the patch boundary. Triangles whose centroids
fall inside the polygon receive `section` instead of the main shell section.

# From Meshes.Ngon (handles any polygon shape, including circular approximations)
```julia
using Meshes
stiff = ShellSection(0.5u"m", 200e9u"Pa", 0.3)
patch = ShellPatch(
    Quadrangle(Point(2.75, 2.75), Point(3.25, 2.75),
               Point(3.25, 3.25), Point(2.75, 3.25)),
    stiff)
```

# From center + dimensions (rectangular shorthand, in meters or Unitful)
```julia
patch = ShellPatch(3.0, 3.0, 0.5, 0.5, stiff)             # meters
patch = ShellPatch(3.0u"m", 3.0u"m", 0.5u"m", 0.5u"m", stiff)  # Unitful
```

# Use with Shell()
```julia
shells = Shell(corners, section; interior_patches=[patch])
```
"""
struct ShellPatch
    vertices::Vector{Tuple{Float64, Float64}}  # polygon vertices (m)
    center::Tuple{Float64, Float64}             # centroid (m)
    section::ShellSection
    id::Symbol

    function ShellPatch(
        vertices::Vector{Tuple{Float64, Float64}},
        center::Tuple{Float64, Float64},
        section::ShellSection,
        id::Symbol = :patch
    )
        @assert length(vertices) >= 3 "ShellPatch requires at least 3 vertices"
        return new(vertices, center, section, id)
    end
end

# Constructor from Meshes.Ngon — auto-computes centroid via Meshes.centroid
function ShellPatch(polygon::Meshes.Ngon, section::ShellSection; id::Symbol = :patch)
    verts = Tuple{Float64, Float64}[]
    for v in Meshes.vertices(polygon)
        c = Meshes.coords(v)
        push!(verts, (_coord_to_meters(c.x), _coord_to_meters(c.y)))
    end
    cen = Meshes.centroid(polygon)
    cc = Meshes.coords(cen)
    center = (_coord_to_meters(cc.x), _coord_to_meters(cc.y))
    return ShellPatch(verts, center, section, id)
end

# Constructor from center + dimensions (rectangular, Float64 in meters)
function ShellPatch(
    cx::Float64, cy::Float64, c1::Float64, c2::Float64,
    section::ShellSection;
    id::Symbol = :patch
)
    verts = [
        (cx - c1/2, cy - c2/2),
        (cx + c1/2, cy - c2/2),
        (cx + c1/2, cy + c2/2),
        (cx - c1/2, cy + c2/2),
    ]
    return ShellPatch(verts, (cx, cy), section, id)
end

# Unitful overload
function ShellPatch(
    cx::Quantity, cy::Quantity, c1::Quantity, c2::Quantity,
    section::ShellSection;
    id::Symbol = :patch
)
    return ShellPatch(
        Float64(ustrip(u"m", cx)), Float64(ustrip(u"m", cy)),
        Float64(ustrip(u"m", c1)), Float64(ustrip(u"m", c2)),
        section; id=id)
end

"""
    _add_patch_geometry!(all_points, patches, boundary_pts) -> Set{Tuple{Int,Int}}

Add patch vertex points to the mesh point list and return constrained segment
index pairs for `DT.triangulate(...; segments=...)`.
"""
function _add_patch_geometry!(
    all_points::Vector{Tuple{Float64, Float64}},
    patches::Vector{ShellPatch},
    boundary_pts::Vector{Tuple{Float64, Float64}}
)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)

    # Build existing-point index for deduplication
    existing = Dict{Tuple{Int64, Int64}, Int}()
    for (i, pt) in enumerate(all_points)
        existing[(round_coord(pt[1]), round_coord(pt[2]))] = i
    end

    segments = Set{Tuple{Int, Int}}()

    for patch in patches
        vertex_indices = Int[]

        for v in patch.vertices
            key = (round_coord(v[1]), round_coord(v[2]))
            if haskey(existing, key)
                push!(vertex_indices, existing[key])
            else
                if _point_inside_polygon(v, boundary_pts) ||
                   _is_on_boundary_edge(v, boundary_pts, tol)
                    push!(all_points, v)
                    idx = length(all_points)
                    existing[key] = idx
                    push!(vertex_indices, idx)
                end
            end
        end

        # Constrained segments between consecutive vertices
        nv = length(vertex_indices)
        for i in 1:nv
            j = mod1(i + 1, nv)
            a, b = vertex_indices[i], vertex_indices[j]
            a != b && push!(segments, (min(a, b), max(a, b)))
        end
    end

    return segments
end

"""
    _apply_patches!(shells, patches)

Reassign section properties for shell elements whose centroids lie inside a patch.
"""
function _apply_patches!(shells::Vector{ShellTri3}, patches::Vector{ShellPatch})
    isempty(patches) && return shells

    for shell in shells
        cx = sum(ustrip(u"m", n.position[1]) for n in shell.nodes) / 3.0
        cy = sum(ustrip(u"m", n.position[2]) for n in shell.nodes) / 3.0

        for patch in patches
            if _point_inside_polygon((cx, cy), patch.vertices)
                shell.thickness = patch.section.thickness
                shell.E = patch.section.E
                shell.ν = patch.section.ν
                shell.ρ = patch.section.ρ
                shell.id = patch.id
                break  # first matching patch wins
            end
        end
    end

    return shells
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
    all_points, support_point_indices, _ = _generate_mesh_points_with_supports(boundary_pts, n, support_lines)
    
    # Build constrained triangulation
    boundary_loop = vcat(collect(1:nc), [1])
    return DT.triangulate(all_points; boundary_nodes=[boundary_loop])
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
    edge_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    interior_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    target_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_radius::Union{Quantity, Real, Nothing} = nothing,
    refinement_targets::Union{Symbol, Vector{Node}, Nothing} = :auto
) where N
    nc = length(corners)
    
    # ===========================================================================
    # Input Validation
    # ===========================================================================
    @assert nc >= 3 "Shell requires at least 3 corner nodes, got $nc"
    
    # Extract corner positions for validation
    boundary_pts = [(ustrip(u"m", node.position[1]), ustrip(u"m", node.position[2])) 
                    for node in corners]
    
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
    
    # Convert interior supports
    support_lines = [_to_support_line(s) for s in interior_supports]
    
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
        _warn_mesh_density(h_far, Lx, Ly)
    end
    
    # ===========================================================================
    # Decide meshing strategy
    # ===========================================================================
    # Structured mesh only for uniform rectangular panels (no refinement).
    # When refinement is requested, Delaunay + rings is always better because it
    # refines radially/locally instead of creating axis-aligned dense bands.
    use_structured = h_far !== nothing && h_near === nothing && _is_rectangular(boundary_pts) && isempty(support_lines) && isempty(interior_nodes) && isempty(interior_patches)
    
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
            
            boundary_loop = vcat(collect(1:nc), [1])
            tri = DT.triangulate(seed_points; boundary_nodes=[boundary_loop])
            
            # Ruppert refinement: spatially-varying area constraint
            far_h  = h_far !== nothing ? h_far : max(Lx, Ly) / max(effective_n, 1)
            A_far  = 0.433 * far_h^2
            A_near = 0.433 * h_near^2
            r_transition = r_trans !== nothing ? r_trans : 5.0 * h_near

            function _needs_refine(_tri, T)
                i, j, k = DT.triangle_vertices(T)
                p, q, r = DT.get_point(_tri, i, j, k)
                A = abs(DT.triangle_area(p, q, r))
                cx = (p[1] + q[1] + r[1]) / 3
                cy = (p[2] + q[2] + r[2]) / 3
                d_min = minimum(hypot(cx - t[1], cy - t[2]) for t in refine_pts)
                α = clamp(d_min / r_transition, 0.0, 1.0)
                A_limit = A_near + α * (A_far - A_near)
                return A > A_limit
            end

            max_pts = max(5000, 20 * length(seed_points))
            tri_backup = deepcopy(tri)
            try
                DT.refine!(tri; min_angle=30.0, custom_constraint=_needs_refine,
                           max_points=max_pts)
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
                    fixity = if i in edge_node_indices
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
        all_points, support_point_indices, interior_node_map = _generate_mesh_points_with_supports(
            boundary_pts, effective_n, support_lines, interior_node_pts)
        
        # Add patch vertices as regular mesh points (no constrained segments).
        # _apply_patches! assigns stiffened sections by centroid check, so
        # strict edge conformity is not required.
        if !isempty(interior_patches)
            _add_patch_geometry!(all_points, interior_patches, boundary_pts)
        end
        
        boundary_loop = vcat(collect(1:nc), [1])
        tri = DT.triangulate(all_points; boundary_nodes=[boundary_loop])
        
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
            
            if i in support_point_indices
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
    edge_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    interior_support_type::Union{Symbol, Vector{Bool}} = :pinned,
    target_edge_length::Union{Quantity, Real, Nothing} = 0.25u"m",
    refinement_edge_length::Union{Quantity, Real, Nothing} = nothing,
    refinement_radius::Union{Quantity, Real, Nothing} = nothing,
    refinement_targets::Union{Symbol, Vector{Node}, Nothing} = :auto
) where N
    return Shell(corners, n, section; 
        id=id, 
        interior_supports=interior_supports,
        interior_nodes=interior_nodes,
        interior_patches=interior_patches,
        edge_support_type=edge_support_type,
        interior_support_type=interior_support_type,
        target_edge_length=target_edge_length,
        refinement_edge_length=refinement_edge_length,
        refinement_radius=refinement_radius,
        refinement_targets=refinement_targets
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
Generate mesh points including support line points and interior nodes.
Returns (all_points, support_point_indices, interior_node_point_indices).

The third return value maps input interior node index → output point index,
allowing the caller to reuse actual Node objects instead of creating new ones.
"""
function _generate_mesh_points_with_supports(
    boundary_pts::Vector{Tuple{Float64, Float64}}, 
    n::Int,
    support_lines::Vector{_SupportLine},
    interior_node_pts::Vector{Tuple{Float64, Float64}} = Tuple{Float64, Float64}[]
)
    nc = length(boundary_pts)
    
    # Tolerance for deduplication (in meters)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)
    
    # Use a dict to dedupe: key is rounded coords, value is (exact coords, is_support, interior_node_idx)
    # interior_node_idx=0 means not an interior node
    points_dict = Dict{Tuple{Int64, Int64}, Tuple{Tuple{Float64, Float64}, Bool, Int}}()
    
    # Add boundary points first (these take priority)
    for pt in boundary_pts
        key = (round_coord(pt[1]), round_coord(pt[2]))
        points_dict[key] = (pt, false, 0)
    end
    
    # Add interior node points (these are actual nodes we want to embed)
    for (idx, pt) in enumerate(interior_node_pts)
        key = (round_coord(pt[1]), round_coord(pt[2]))
        if !haskey(points_dict, key)
            # Only add if inside polygon or on boundary
            if _point_inside_polygon(pt, boundary_pts) || _is_on_boundary_edge(pt, boundary_pts, tol)
                points_dict[key] = (pt, false, idx)  # Track which interior node
            end
        end
    end
    
    # Add support line points
    for line in support_lines
        line_pts = _points_along_line(line.start_pt, line.end_pt, n)
        for pt in line_pts
            key = (round_coord(pt[1]), round_coord(pt[2]))
            if !haskey(points_dict, key)
                # Only add if inside polygon
                if _point_inside_polygon(pt, boundary_pts) || _is_on_boundary_edge(pt, boundary_pts, tol)
                    points_dict[key] = (pt, true, 0)  # Mark as support point
                end
            end
        end
    end
    
    if n <= 1
        # Build result with support tracking
        result = collect(boundary_pts)
        support_indices = Set{Int}()
        interior_node_map = Dict{Int, Int}()  # input idx → point idx
        
        for (key, (pt, is_support, interior_idx)) in points_dict
            is_boundary = any(abs(pt[1] - bp[1]) < tol && abs(pt[2] - bp[2]) < tol for bp in boundary_pts)
            if !is_boundary
                push!(result, pt)
                point_idx = length(result)
                if is_support
                    push!(support_indices, point_idx)
                end
                if interior_idx > 0
                    interior_node_map[interior_idx] = point_idx
                end
            end
        end
        
        return result, support_indices, interior_node_map
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
                points_dict[key] = (pt, false, 0)
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
                    points_dict[key] = (pt, false, 0)
                end
            end
        end
    end
    
    # Build result: boundary points first, then others
    result = collect(boundary_pts)
    support_indices = Set{Int}()
    interior_node_map = Dict{Int, Int}()  # input idx → point idx
    
    for (key, (pt, is_support, interior_idx)) in points_dict
        is_boundary = any(abs(pt[1] - bp[1]) < tol && abs(pt[2] - bp[2]) < tol for bp in boundary_pts)
        if !is_boundary
            push!(result, pt)
            point_idx = length(result)
            if is_support
                push!(support_indices, point_idx)
            end
            if interior_idx > 0
                interior_node_map[interior_idx] = point_idx
            end
        end
    end
    
    return result, support_indices, interior_node_map
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

# ============================================================================
# Structured mesh → ShellTri3 elements
# ============================================================================

"""
    _create_shells_from_structured(xs, ys, corners, boundary_pts, interior_nodes,
                                    section, id, edge_support_type, interior_support_type)

Build ShellTri3 elements from a structured rectangular grid via `_t3blockx`.
Maps grid nodes to Asap `Node`s with appropriate fixities.
"""
function _create_shells_from_structured(
    xs::Vector{Float64},
    ys::Vector{Float64},
    corners,
    boundary_pts::Vector{Tuple{Float64, Float64}},
    interior_nodes::Vector{Node},
    section::ShellSection,
    id::Symbol,
    edge_support_type::Union{Symbol, Vector{Bool}},
    interior_support_type::Union{Symbol, Vector{Bool}}
)
    # Generate structured mesh
    mesh_pts, mesh_conn = _t3blockx(xs, ys)
    
    # Average z from corners
    nc = length(corners)
    z_val = sum(ustrip(u"m", node.position[3]) for node in corners) / nc
    
    tol = 1e-6
    
    # Build map from interior node positions to actual Node objects (for structural connectivity)
    interior_map = Dict{Tuple{Int64, Int64}, Node}()
    round_coord(x) = round(Int64, x / tol)
    for nd in interior_nodes
        key = (round_coord(ustrip(u"m", nd.position[1])),
               round_coord(ustrip(u"m", nd.position[2])))
        interior_map[key] = nd
    end
    
    # Build map from corner positions to actual corner Node objects
    corner_map = Dict{Tuple{Int64, Int64}, Node}()
    for corner in corners
        key = (round_coord(ustrip(u"m", corner.position[1])),
               round_coord(ustrip(u"m", corner.position[2])))
        corner_map[key] = corner
    end
    
    # Create nodes for each mesh point
    node_array = Vector{Node}(undef, length(mesh_pts))
    
    for (i, pt) in enumerate(mesh_pts)
        key = (round_coord(pt[1]), round_coord(pt[2]))
        
        if haskey(corner_map, key)
            # Use the actual corner node
            node_array[i] = corner_map[key]
        elseif haskey(interior_map, key)
            # Use the actual interior node (structural connectivity)
            node_array[i] = interior_map[key]
        else
            pos = [pt[1], pt[2], z_val] .* u"m"
            
            # Determine fixity
            is_edge = _is_on_boundary_edge(pt, boundary_pts, tol)
            fixity = is_edge ? edge_support_type : :free
            
            node_array[i] = _create_node_with_fixity(pos, fixity, id)
        end
    end
    
    # Create shell elements
    thickness_q = section.thickness * u"m"
    shells = ShellTri3[]
    for (i, j, k) in mesh_conn
        push!(shells, ShellTri3(
            (node_array[i], node_array[j], node_array[k]),
            thickness_q,
            section.E * u"Pa",
            section.ν;
            ρ=section.ρ,
            id=id
        ))
    end
    
    return shells
end

# ============================================================================
# Refinement rings for Delaunay meshing
# ============================================================================

"""
    _add_refinement_rings!(all_points, refine_pts, h_near, r_transition,
                           boundary_pts, patches=ShellPatch[])

Add concentric rings of points around each refinement target to improve
Delaunay mesh density near supports.

Points are discarded if they:
- fall outside the slab boundary polygon,
- duplicate an existing mesh point, or
- fall **inside any ShellPatch** interior (column footprint).  Placing ring
  points inside a patch creates degenerate slivers between the ring nodes
  and the patch's constrained edges, tanking the stiffness-matrix condition
  number.  Patch geometry must be added to `all_points` *before* calling
  this function.
"""
function _add_refinement_rings!(
    all_points::Vector{Tuple{Float64, Float64}},
    refine_pts::Vector{Tuple{Float64, Float64}},
    h_near::Float64,
    r_transition::Float64,
    boundary_pts::Vector{Tuple{Float64, Float64}},
    patches::Vector{ShellPatch} = ShellPatch[]
)
    tol = 1e-6
    round_coord(x) = round(Int64, x / tol)
    
    # Existing point set for deduplication
    existing = Set{Tuple{Int64, Int64}}()
    for pt in all_points
        push!(existing, (round_coord(pt[1]), round_coord(pt[2])))
    end

    # Precompute patch vertex arrays for point-in-polygon tests
    patch_polys = [p.vertices for p in patches]
    
    for center in refine_pts
        # Add rings at increasing radii
        radii = Float64[]
        r = h_near
        while r < r_transition
            push!(radii, r)
            r *= 1.6  # growth factor
        end
        
        for r in radii
            # Number of points on this ring (circumference / desired spacing)
            n_ring = max(6, round(Int, 2π * r / h_near))
            for k in 0:n_ring-1
                θ = 2π * k / n_ring
                pt = (center[1] + r * cos(θ), center[2] + r * sin(θ))
                key = (round_coord(pt[1]), round_coord(pt[2]))
                
                # Skip duplicates
                key in existing && continue
                # Must be inside the slab boundary
                _point_inside_polygon(pt, boundary_pts) || continue
                # Must NOT be inside any patch interior
                any(poly -> _point_inside_polygon(pt, poly), patch_polys) && continue
                
                push!(all_points, pt)
                push!(existing, key)
            end
        end
    end
end
