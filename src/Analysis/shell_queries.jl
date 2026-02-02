#=
Shell Spatial Queries and Region Integration
=============================================

Functions for querying shell elements by location and integrating
results over polygonal regions.

Key functions:
- `shell_tris_at_point`: Find shell triangles at/near a point
- `shell_tris_in_region`: Find shell triangles within a polygon
- `bending_moments(model; polygon=...)`: Integrate moments over a region
- `bending_moments(model; pts=...)`: Query moments at multiple points
=#

# ============================================================================
# Geometry Helpers
# ============================================================================

"""
    shell_centroid(elem::ShellElement) -> NamedTuple{(:x, :y)}

Return (x=..., y=...) centroid of a shell element in global coordinates (meters).
"""
function shell_centroid(elem::ShellElement)
    x = sum(ustrip(u"m", n.position[1]) for n in elem.nodes) / length(elem.nodes)
    y = sum(ustrip(u"m", n.position[2]) for n in elem.nodes) / length(elem.nodes)
    return (x=x, y=y)
end

"""
    shell_centroid_3d(elem::ShellElement) -> NamedTuple{(:x, :y, :z)}

Return (x=..., y=..., z=...) centroid of a shell element in global coordinates (meters).
"""
function shell_centroid_3d(elem::ShellElement)
    n = length(elem.nodes)
    x = sum(ustrip(u"m", node.position[1]) for node in elem.nodes) / n
    y = sum(ustrip(u"m", node.position[2]) for node in elem.nodes) / n
    z = sum(ustrip(u"m", node.position[3]) for node in elem.nodes) / n
    return (x=x, y=y, z=z)
end

"""
    point_in_triangle_with_tol(pt, elem::ShellTri3, tol) -> Bool

Check if a 2D point is inside or near a triangular shell element.
Uses barycentric coordinates with a tolerance band for edge cases.

# Arguments
- `pt`: Query point as `(x, y)` tuple (meters)
- `elem`: Shell triangle element
- `tol`: Distance tolerance in meters (allows points slightly outside)
"""
function point_in_triangle_with_tol(pt::Tuple{T,T}, elem::ShellTri3, tol::Float64) where T<:Real
    # Get triangle vertices (2D projection)
    v1 = (ustrip(u"m", elem.nodes[1].position[1]), ustrip(u"m", elem.nodes[1].position[2]))
    v2 = (ustrip(u"m", elem.nodes[2].position[1]), ustrip(u"m", elem.nodes[2].position[2]))
    v3 = (ustrip(u"m", elem.nodes[3].position[1]), ustrip(u"m", elem.nodes[3].position[2]))
    
    # Compute barycentric coordinates
    denom = (v2[2] - v3[2]) * (v1[1] - v3[1]) + (v3[1] - v2[1]) * (v1[2] - v3[2])
    
    if abs(denom) < 1e-12
        return false  # Degenerate triangle
    end
    
    λ1 = ((v2[2] - v3[2]) * (pt[1] - v3[1]) + (v3[1] - v2[1]) * (pt[2] - v3[2])) / denom
    λ2 = ((v3[2] - v1[2]) * (pt[1] - v3[1]) + (v1[1] - v3[1]) * (pt[2] - v3[2])) / denom
    λ3 = 1.0 - λ1 - λ2
    
    # Tolerance band: allow slightly outside (negative λ down to -tol/h)
    # where h is characteristic element size
    h = sqrt(elem.area)
    λ_tol = h > 0 ? tol / h : 0.0
    
    return (λ1 >= -λ_tol) && (λ2 >= -λ_tol) && (λ3 >= -λ_tol)
end

"""
    _point_in_polygon(pt, polygon) -> Bool

Check if a 2D point is inside a polygon using ray casting algorithm.

# Arguments
- `pt`: Query point as `(x, y)` tuple or NamedTuple with x, y fields
- `polygon`: Vector of `(x, y)` tuples defining polygon vertices (closed loop)
"""
function _point_in_polygon(pt::Tuple{T,T}, polygon::Vector{Tuple{S,S}}) where {T<:Real, S<:Real}
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

# NamedTuple overload for shell_centroid compatibility
function _point_in_polygon(pt::NamedTuple{(:x, :y)}, polygon::Vector{Tuple{S,S}}) where S<:Real
    return _point_in_polygon((pt.x, pt.y), polygon)
end

# ============================================================================
# Spatial Query Functions
# ============================================================================

"""
    shell_tris_at_point(shells, pt; tol=1e-6) -> Vector{ShellTri3}
    shell_tris_at_point(model, pt; tol=1e-6) -> Vector{ShellTri3}

Find all shell triangles that contain or nearly contain a point.

At element edges/corners, multiple triangles may be returned.
The tolerance allows for floating-point imprecision in queries.

# Arguments
- `shells`: Vector of shell elements to search
- `model`: Model containing shell elements (uses `model.shell_elements`)
- `pt`: Query point as `(x, y)` tuple (meters)
- `tol`: Distance tolerance in meters (default: 1e-6 = 1 micron)

# Returns
Vector of `ShellTri3` elements (may be empty, one, or multiple).

# Example
```julia
# Find shells at a column location
shells = shell_tris_at_point(model, (5.0, 3.0); tol=0.01)
```
"""
function shell_tris_at_point(shells::Vector{<:ShellElement}, pt::Tuple{T,T}; 
                              tol::Float64=1e-6) where T<:Real
    result = ShellTri3[]
    
    for shell in shells
        if shell isa ShellTri3 && point_in_triangle_with_tol(pt, shell, tol)
            push!(result, shell)
        end
    end
    
    return result
end

function shell_tris_at_point(model::AbstractModel, pt::Tuple{T,T}; tol::Float64=1e-6) where T<:Real
    shells = _get_shell_elements(model)
    return shell_tris_at_point(shells, pt; tol=tol)
end

# Helper to get shell elements from any model type
_get_shell_elements(model::Model) = model.shell_elements
_get_shell_elements(model::ElementModel) = model.elements  # ShellModel, etc.

"""
    shell_tris_in_region(shells, polygon) -> Vector{ShellTri3}
    shell_tris_in_region(model, polygon) -> Vector{ShellTri3}

Find all shell triangles whose centroids fall within a polygon boundary.

# Arguments
- `shells`: Vector of shell elements to search
- `model`: Model containing shell elements
- `polygon`: Vector of `(x, y)` tuples defining the region boundary

# Returns
Vector of `ShellTri3` elements within the region.

# Example
```julia
# Get shells within a column strip
strip_boundary = [(0.0, 0.0), (10.0, 0.0), (10.0, 2.5), (0.0, 2.5)]
tris = shell_tris_in_region(model, strip_boundary)
```
"""
function shell_tris_in_region(shells::Vector{<:ShellElement}, 
                               polygon::Vector{Tuple{T,T}}) where T<:Real
    result = ShellTri3[]
    
    for shell in shells
        if shell isa ShellTri3
            centroid = shell_centroid(shell)
            if _point_in_polygon(centroid, polygon)
                push!(result, shell)
            end
        end
    end
    
    return result
end

function shell_tris_in_region(model::AbstractModel, polygon::Vector{Tuple{T,T}}) where T<:Real
    shells = _get_shell_elements(model)
    return shell_tris_in_region(shells, polygon)
end

# ============================================================================
# Bending Moments - Extended Dispatch
# ============================================================================

# Note: The basic bending_moments(elem, u) is defined in shell.jl
# Here we add dispatch for:
# - Single element + model (convenience)
# - Vector of elements + model
# - Single point query
# - Polygon region integration
# - Multiple point queries
# - Element subset + region/point queries

"""
    bending_moments(elem::ShellTri3, model::AbstractModel) -> Vector{Float64}

Compute bending moments [Mxx, Myy, Mxy] for a single element using model's displacement.
"""
function bending_moments(elem::ShellTri3, model::AbstractModel)
    @assert model.processed "Model must be solved first (call solve!)"
    return bending_moments(elem, model.u)
end

"""
    bending_moments(model::Model, pt::Tuple{T,T}; tol=1e-6) -> Union{Vector{Float64}, Nothing}

Query bending moments at a single point.
If multiple triangles cover the point (edge/corner), returns averaged moments.

# Arguments
- `model`: Solved model
- `pt`: Query point as `(x, y)` tuple (meters)
- `tol`: Tolerance for point-in-triangle test (meters)

# Returns
`[Mxx, Myy, Mxy]` or `nothing` if no shells found at point.
"""
function bending_moments(model::AbstractModel, pt::Tuple{T,T}; tol::Float64=1e-6) where T<:Real
    @assert model.processed "Model must be solved first (call solve!)"
    
    tris = shell_tris_at_point(model, pt; tol=tol)
    
    if isempty(tris)
        return nothing
    elseif length(tris) == 1
        return bending_moments(tris[1], model.u)
    else
        # Average moments from all touching elements
        moments = [bending_moments(t, model.u) for t in tris]
        return [sum(m[i] for m in moments) / length(moments) for i in 1:3]
    end
end

"""
    bending_moments(model::Model; polygon=nothing, pts=nothing, tol=1e-6)

Query bending moments by region (polygon) or at multiple points (pts).

# Keyword Arguments
- `polygon`: Vector of `(x, y)` tuples defining a region for integration
- `pts`: Vector of `(x, y)` tuples for point queries
- `tol`: Tolerance for point queries (meters, default: 1e-6)

Only one of `polygon` or `pts` should be specified.

# Returns

For `polygon`: NamedTuple with:
- `Mxx, Myy, Mxy`: Area-weighted total moment (N·m)
- `Mxx_avg, Myy_avg, Mxy_avg`: Average moment intensity (N·m/m)
- `Mxx_max, Myy_max, Mxy_max`: Peak absolute values
- `area`: Total shell area in region (m²)
- `shell_tris`: Vector of shells in the region

For `pts`: Vector of `[Mxx, Myy, Mxy]` (or `nothing` for points with no shells)

# Example
```julia
# Region integration
strip = [(0.0, 0.0), (10.0, 0.0), (10.0, 2.5), (0.0, 2.5)]
result = bending_moments(model; polygon=strip)
println("Peak Mxx = \$(result.Mxx_max) N·m/m")

# Point queries
pts = [(1.0, 1.0), (5.0, 5.0), (9.0, 9.0)]
moments = bending_moments(model; pts=pts, tol=0.01)
```
"""
function bending_moments(model::AbstractModel; 
                         polygon::Union{Nothing, Vector{Tuple{T,T}}}=nothing,
                         pts::Union{Nothing, Vector{Tuple{T,T}}}=nothing,
                         tol::Float64=1e-6) where T<:Real
    @assert model.processed "Model must be solved first (call solve!)"
    @assert !isnothing(polygon) ⊻ !isnothing(pts) "Specify exactly one of polygon= or pts="
    
    shells = _get_shell_elements(model)
    if !isnothing(polygon)
        return _integrate_bending_moments(shells, model.u, polygon)
    else
        return [bending_moments(model, pt; tol=tol) for pt in pts]
    end
end

"""
    bending_moments(elems::Vector{<:ShellElement}, model::Model; polygon=nothing, pts=nothing, tol=1e-6)

Query bending moments from a specific set of elements.

# Usage

Without keyword arguments: Returns moments for each element in the vector.
With `polygon=`: Integrates over shells whose centroids fall within the polygon.
With `pts=`: Queries at each point, averaging over touching elements.

# Example
```julia
# Get moments for all elements in a vector
Ms = bending_moments(shells[1:5], model)

# Get shells for a specific slab, then integrate over a strip
slab_tris = shell_tris_in_region(model, slab_boundary)
result = bending_moments(slab_tris, model; polygon=strip_boundary)
```
"""
function bending_moments(elems::Vector{<:ShellElement}, model::AbstractModel;
                         polygon::Union{Nothing, Vector{Tuple{T,T}}}=nothing,
                         pts::Union{Nothing, Vector{Tuple{T,T}}}=nothing,
                         tol::Float64=1e-6) where T<:Real
    @assert model.processed "Model must be solved first (call solve!)"
    
    # If neither polygon nor pts is specified, return moments for all elements
    if isnothing(polygon) && isnothing(pts)
        return [bending_moments(e, model.u) for e in elems if e isa ShellTri3]
    end
    
    @assert !isnothing(polygon) ⊻ !isnothing(pts) "Specify at most one of polygon= or pts="
    
    if !isnothing(polygon)
        return _integrate_bending_moments(elems, model.u, polygon)
    else
        return [_query_moments_at_point(elems, model.u, pt, tol) for pt in pts]
    end
end

# ============================================================================
# Internal Integration/Query Functions
# ============================================================================

"""
Integrate bending moments over shells within a polygon region.
"""
function _integrate_bending_moments(shells::Vector{<:ShellElement}, 
                                     u::Vector{Float64},
                                     polygon::Vector{Tuple{T,T}}) where T<:Real
    tris = shell_tris_in_region(shells, polygon)
    
    if isempty(tris)
        return (
            Mxx = 0.0, Myy = 0.0, Mxy = 0.0,
            Mxx_avg = 0.0, Myy_avg = 0.0, Mxy_avg = 0.0,
            Mxx_max = 0.0, Myy_max = 0.0, Mxy_max = 0.0,
            area = 0.0,
            shell_tris = tris
        )
    end
    
    Mxx_total, Myy_total, Mxy_total = 0.0, 0.0, 0.0
    Mxx_max, Myy_max, Mxy_max = -Inf, -Inf, -Inf
    area_total = 0.0
    
    for tri in tris
        M = bending_moments(tri, u)  # [Mxx, Myy, Mxy] in N·m/m
        A = tri.area  # m²
        
        # Weighted accumulation (M is moment per unit width, M*A = total moment over element)
        Mxx_total += M[1] * A
        Myy_total += M[2] * A
        Mxy_total += M[3] * A
        area_total += A
        
        # Track peaks (absolute values for max)
        Mxx_max = max(Mxx_max, abs(M[1]))
        Myy_max = max(Myy_max, abs(M[2]))
        Mxy_max = max(Mxy_max, abs(M[3]))
    end
    
    return (
        Mxx = Mxx_total,
        Myy = Myy_total,
        Mxy = Mxy_total,
        Mxx_avg = Mxx_total / area_total,
        Myy_avg = Myy_total / area_total,
        Mxy_avg = Mxy_total / area_total,
        Mxx_max = Mxx_max,
        Myy_max = Myy_max,
        Mxy_max = Mxy_max,
        area = area_total,
        shell_tris = tris
    )
end

"""
Query moments at a point from a specific element set.
"""
function _query_moments_at_point(shells::Vector{<:ShellElement}, u::Vector{Float64},
                                  pt::Tuple{T,T}, tol::Float64) where T<:Real
    tris = shell_tris_at_point(shells, pt; tol=tol)
    
    if isempty(tris)
        return nothing
    elseif length(tris) == 1
        return bending_moments(tris[1], u)
    else
        moments = [bending_moments(t, u) for t in tris]
        return [sum(m[i] for m in moments) / length(moments) for i in 1:3]
    end
end
