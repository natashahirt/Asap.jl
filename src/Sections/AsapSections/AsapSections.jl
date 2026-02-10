module AsapSections

using LinearAlgebra

# ── Polygonal section types ──────────────────────────────────────────
include("Sections.jl")
export AbstractPolygonalSection, PolygonalSection
export SolidSection, VoidSection, CompoundSection, OffsetSection

# ── Extended section properties (Zx, Zy, rx, ry) ────────────────────
include("SectionProperties.jl")
export SectionProperties

# ── Geometry utilities ───────────────────────────────────────────────
include("GeneralFunctions.jl")
export poly_area
export center_at_centroid!
export rotate_section!
export translate_section!

# ── Depth / compression-zone analysis (Sutherland-Hodgman) ──────────
include("DepthAnalysis.jl")
export sutherland_hodgman
export sutherland_hodgman_abs
export intersection
export depth_map, depth_map_abs
export area_from_depth, area_from_depth_abs
export depth_from_area

end # module AsapSections
