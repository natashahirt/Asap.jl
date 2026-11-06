
# =============================================================================
# Internal Helpers
# =============================================================================

"""
Internal: Generate TributaryLoads from multiple shell panels.

For multi-bay slabs where each bay is a separate panel (Vector of shells),
this function combines tributary contributions from adjacent panels.
Interior beams that appear in multiple panels get combined load.
"""
function _shell_panels_to_tributary_loads(
    panels::Vector{Vector{S}},
    beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
    tol::Float64 = 1e-6
)::Vector{TributaryLoad} where S <: ShellElement
    isempty(panels) && return TributaryLoad[]
    isempty(beams) && return TributaryLoad[]
    
    # Track which beams are matched by at least one panel
    beam_matched = Dict{FrameElement, Bool}(b => false for b in beams)
    
    # Accumulate tributary data per beam: beam → list of (positions, widths) tuples
    beam_tributaries = Dict{FrameElement, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}()
    
    for panel_shells in panels
        isempty(panel_shells) && continue
        
        # Extract panel boundary
        boundary = _get_shell_panel_boundary(panel_shells)
        isempty(boundary) && continue
        
        boundary_points = [Point(p[1]*u"m", p[2]*u"m") for p in boundary]
        
        # Compute tributary polygons for this panel
        trib_polys = isnothing(axis) ? 
            get_tributary_polygons(boundary_points) :
            get_tributary_polygons(boundary_points; axis=axis)
        
        # Match beams to this panel's edges
        for beam in beams
            edge_idx = _match_beam_to_polygon_edge(beam, boundary_points; tol=tol)
            if !isnothing(edge_idx)
                beam_matched[beam] = true
                tp = trib_polys[edge_idx]
                
                # Extract depth profile from polygon vertices
                positions, depths = _extract_depth_profile(tp)
                
                if !haskey(beam_tributaries, beam)
                    beam_tributaries[beam] = Tuple{Vector{Float64}, Vector{Float64}}[]
                end
                push!(beam_tributaries[beam], (positions, depths))
            end
        end
    end
    
    # Check for unmatched beams
    unmatched = [b for (b, matched) in beam_matched if !matched]
    if !isempty(unmatched)
        beam_ids = [b.id for b in unmatched]
        error("Beams do not match any panel edge: $beam_ids. " *
              "All beams must lie on at least one panel boundary.")
    end
    
    # Create combined TributaryLoads
    loads = TributaryLoad[]
    
    for beam in beams
        tributaries = get(beam_tributaries, beam, Tuple{Vector{Float64}, Vector{Float64}}[])
        
        if isempty(tributaries)
            continue
        elseif length(tributaries) == 1
            positions, depths = tributaries[1]
        else
            # Multiple panels - combine by summing depths
            combined = _combine_tributary_depths(tributaries)
            positions, depths = combined
        end
        
        # Clamp negative depths to zero
        depths = max.(depths, 0.0)
        widths = [d * u"m" for d in depths]
        
        # Handle zero-load case
        max_depth = isempty(depths) ? 0.0 : maximum(depths)
        if max_depth < 1e-9
            positions = [0.0, 1.0]
            widths = [0.0u"m", 0.0u"m"]
        end
        
        push!(loads, TributaryLoad(beam, positions, widths, pressure, direction))
    end
    
    return loads
end

"""Internal: Combine multiple tributary depth profiles by summing depths at each position."""
function _combine_tributary_depths(
    tributaries::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
)::Tuple{Vector{Float64}, Vector{Float64}}
    # Collect all unique positions
    all_positions = Set{Float64}()
    for (positions, _) in tributaries
        union!(all_positions, positions)
    end
    
    sorted_positions = sort(collect(all_positions))
    
    # Sum depths at each position (interpolate if needed)
    combined_depths = zeros(length(sorted_positions))
    
    for (positions, depths) in tributaries
        for (i, s) in enumerate(sorted_positions)
            combined_depths[i] += _interpolate_depth(s, positions, depths)
        end
    end
    
    return (sorted_positions, combined_depths)
end

"""
Linear interpolation of depth at position s.
"""
function _interpolate_depth(s::Float64, positions::Vector{Float64}, depths::Vector{Float64})::Float64
    n = length(positions)
    n == 0 && return 0.0
    n == 1 && return depths[1]
    
    # Clamp to range
    s <= positions[1] && return depths[1]
    s >= positions[end] && return depths[end]
    
    # Find bracketing interval
    for i in 1:(n-1)
        if positions[i] <= s <= positions[i+1]
            t = (s - positions[i]) / (positions[i+1] - positions[i])
            return depths[i] + t * (depths[i+1] - depths[i])
        end
    end
    
    return 0.0
end
