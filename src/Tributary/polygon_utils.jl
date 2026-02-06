# =============================================================================
# Polygon Utility Functions
# =============================================================================
#
# General-purpose polygon geometry utilities.
#
# =============================================================================

import Meshes
using Unitful: ustrip

"""
    is_convex_polygon(vertices) -> Bool

Check if a polygon (given as ordered vertices) is convex.

Uses the cross-product test: all cross products of consecutive edge pairs
must have the same sign for a convex polygon.

# Arguments
- `vertices`: Vector of coordinates. Accepts multiple formats:
  - Tuples: `[(0.0, 0.0), (1.0, 0.0), ...]`
  - Vectors: `[[0.0, 0.0], [1.0, 0.0], ...]`
  - Meshes.Point: `[Point(0.0, 0.0), ...]`
  - Any type with `.x` and `.y` fields or `[1]`/`[2]` indexing

# Returns
- `true` if polygon is convex, `false` if concave

# Examples
```julia
# Square (convex)
is_convex_polygon([(0, 0), (1, 0), (1, 1), (0, 1)])  # true

# L-shape (concave)
is_convex_polygon([(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)])  # false
```
"""
function is_convex_polygon(vertices::Vector)
    n = length(vertices)
    n < 3 && return true  # Degenerate cases are "convex"
    
    sign_pos = 0
    sign_neg = 0
    
    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[mod1(i + 1, n)]
        p3 = vertices[mod1(i + 2, n)]
        
        # Edge vectors
        dx1 = _poly_get_x(p2) - _poly_get_x(p1)
        dy1 = _poly_get_y(p2) - _poly_get_y(p1)
        dx2 = _poly_get_x(p3) - _poly_get_x(p2)
        dy2 = _poly_get_y(p3) - _poly_get_y(p2)
        
        # Cross product (z-component)
        cross = dx1 * dy2 - dy1 * dx2
        
        if cross > 1e-10
            sign_pos += 1
        elseif cross < -1e-10
            sign_neg += 1
        end
        # Collinear points (cross ≈ 0) don't affect convexity
    end
    
    # Convex if all cross products have same sign (or are zero)
    return sign_pos == 0 || sign_neg == 0
end

# =============================================================================
# Coordinate Accessors (support multiple input formats)
# Prefixed with _poly_ to avoid name collisions with other Asap modules
# =============================================================================

# Tuples
_poly_get_x(p::Tuple) = Float64(p[1])
_poly_get_y(p::Tuple) = Float64(p[2])

# AbstractVector (includes StaticArrays.SVector)
_poly_get_x(p::AbstractVector) = Float64(p[1])
_poly_get_y(p::AbstractVector) = Float64(p[2])

# Named tuple or struct with .x, .y (fallback)
_poly_get_x(p) = Float64(p.x)
_poly_get_y(p) = Float64(p.y)

# Meshes.Point - extract via coords
function _poly_get_x(p::Meshes.Point)
    c = Meshes.coords(p)
    Float64(ustrip(c.x))
end

function _poly_get_y(p::Meshes.Point)
    c = Meshes.coords(p)
    Float64(ustrip(c.y))
end
