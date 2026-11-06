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
    _add_constrained_line_geometry!(all_points, lines, boundary_pts) -> Set{Tuple{Int,Int}}

Add interior constrained-line polyline vertices to the mesh point list and
return constrained segment index pairs for `DT.triangulate(...; segments=...)`.

Each entry of `lines` is an **open** polyline `[(x1,y1), (x2,y2), ...]`; segments
connect consecutive vertices only (no wrap-around edge, unlike `_add_patch_geometry!`
which closes the polygon). The returned tuples are `(min(a,b), max(a,b))` index
pairs indexing into `all_points`.

Vertices outside the boundary polygon (and not on its edge) are dropped; the
polyline is simply broken there (the two segments touching a dropped vertex are
omitted). Callers must ensure no other mesh point lies on a segment interior, or
`DT.triangulate` will split the constrained edge at that point.
"""
function _add_constrained_line_geometry!(
    all_points::Vector{Tuple{Float64, Float64}},
    lines::Vector{Vector{Tuple{Float64, Float64}}},
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

    for line in lines
        # Resolve each vertex to a point index (add if missing & inside/on boundary).
        # Index 0 is a sentinel for a vertex dropped as outside the domain.
        vertex_indices = Int[]
        for v in line
            key = (round_coord(v[1]), round_coord(v[2]))
            if haskey(existing, key)
                push!(vertex_indices, existing[key])
            elseif _point_inside_polygon(v, boundary_pts) ||
                   _is_on_boundary_edge(v, boundary_pts, tol)
                push!(all_points, v)
                idx = length(all_points)
                existing[key] = idx
                push!(vertex_indices, idx)
            else
                push!(vertex_indices, 0)
            end
        end

        # Constrained segments between consecutive in-domain vertices
        for i in 1:(length(vertex_indices) - 1)
            a, b = vertex_indices[i], vertex_indices[i + 1]
            (a == 0 || b == 0 || a == b) && continue
            push!(segments, (min(a, b), max(a, b)))
        end
    end

    return segments
end

"""
    _on_any_constrained_segment_interior(pt, lines, tol) -> Bool

True if `pt` lies on the interior (not within `tol` of a vertex) of any polyline
segment in `lines`. Used to drop auto-generated mesh points that would otherwise
force `DT.triangulate` to split a constrained edge.
"""
function _on_any_constrained_segment_interior(
    pt::Tuple{Float64, Float64},
    lines::Vector{Vector{Tuple{Float64, Float64}}},
    tol::Float64
)
    for line in lines
        for i in 1:(length(line) - 1)
            a, b = line[i], line[i + 1]
            hypot(pt[1] - a[1], pt[2] - a[2]) < tol && continue
            hypot(pt[1] - b[1], pt[2] - b[2]) < tol && continue
            _point_on_segment(pt, a, b, tol) && return true
        end
    end
    return false
end

"""
    _sample_constrained_line_points(lines, spacing) -> Vector{Tuple{Float64,Float64}}

Sample points along each polyline segment at roughly `spacing` intervals
(endpoints included). Used as refinement targets so triangle quality stays good
near constrained edges when Ruppert refinement is active.
"""
function _sample_constrained_line_points(
    lines::Vector{Vector{Tuple{Float64, Float64}}},
    spacing::Float64
)
    pts = Tuple{Float64, Float64}[]
    spacing <= 0 && return pts
    for line in lines
        for i in 1:(length(line) - 1)
            a, b = line[i], line[i + 1]
            L = hypot(b[1] - a[1], b[2] - a[2])
            L < 1e-12 && continue
            nseg = max(1, ceil(Int, L / spacing))
            for k in 0:nseg
                t = k / nseg
                push!(pts, (a[1] * (1 - t) + b[1] * t, a[2] * (1 - t) + b[2] * t))
            end
        end
    end
    return pts
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

