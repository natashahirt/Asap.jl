#=
Shell-to-Beam Load Distribution (Internal Implementation)
=========================================================

Internal functions for converting shell panel geometry into TributaryLoads.
These are implementation details used by AreaLoad when distribute_to=elements.

Public API: Use `AreaLoad(shells, pressure; distribute_to=beams)` instead.

Attribution: Part of the ASAP structural analysis framework.
=#

using Meshes: Point, coords
using Unitful: ustrip, @u_str, Quantity

# =============================================================================
# Edge Matching (Internal)
# =============================================================================

"""Internal: Match beam to polygon edge index."""
function _match_beam_to_polygon_edge(
    beam::FrameElement,
    vertices::Vector{<:Point};
    tol::Float64 = 1e-6
)::Union{Int, Nothing}
    n = length(vertices)
    n < 2 && return nothing
    
    # Get beam endpoint positions in meters (Node stores position as a vector)
    # Element uses nodeStart and nodeEnd fields
    beam_start = (
        Float64(ustrip(u"m", beam.nodeStart.position[1])),
        Float64(ustrip(u"m", beam.nodeStart.position[2]))
    )
    beam_end = (
        Float64(ustrip(u"m", beam.nodeEnd.position[1])),
        Float64(ustrip(u"m", beam.nodeEnd.position[2]))
    )
    
    # Check each polygon edge
    for i in 1:n
        j = mod1(i + 1, n)
        
        # Get edge endpoint positions
        c_i = coords(vertices[i])
        c_j = coords(vertices[j])
        v1 = (Float64(ustrip(u"m", c_i.x)), Float64(ustrip(u"m", c_i.y)))
        v2 = (Float64(ustrip(u"m", c_j.x)), Float64(ustrip(u"m", c_j.y)))
        
        # Check both orientations
        dist_forward = hypot(beam_start[1] - v1[1], beam_start[2] - v1[2]) +
                       hypot(beam_end[1] - v2[1], beam_end[2] - v2[2])
        dist_reverse = hypot(beam_start[1] - v2[1], beam_start[2] - v2[2]) +
                       hypot(beam_end[1] - v1[1], beam_end[2] - v1[2])
        
        if dist_forward < 2 * tol || dist_reverse < 2 * tol
            return i
        end
    end
    
    return nothing
end

"""Internal: Match beam (as tuple coords) to polygon edge index."""
function _match_beam_to_polygon_edge(
    beam_start::NTuple{2, Float64},
    beam_end::NTuple{2, Float64},
    vertices::Vector{NTuple{2, Float64}};
    tol::Float64 = 1e-6
)::Union{Int, Nothing}
    n = length(vertices)
    n < 2 && return nothing
    
    for i in 1:n
        j = mod1(i + 1, n)
        v1, v2 = vertices[i], vertices[j]
        
        dist_forward = hypot(beam_start[1] - v1[1], beam_start[2] - v1[2]) +
                       hypot(beam_end[1] - v2[1], beam_end[2] - v2[2])
        dist_reverse = hypot(beam_start[1] - v2[1], beam_start[2] - v2[2]) +
                       hypot(beam_end[1] - v1[1], beam_end[2] - v1[2])
        
        if dist_forward < 2 * tol || dist_reverse < 2 * tol
            return i
        end
    end
    
    return nothing
end
