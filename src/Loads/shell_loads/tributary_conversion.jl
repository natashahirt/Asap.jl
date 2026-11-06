
# =============================================================================
# Shell to Tributary Load Conversion (Internal)
# =============================================================================

"""Internal: Main implementation for shell-to-tributary conversion."""
function _shell_to_tributary_loads(
    shells::Vector{<:ShellElement},
    beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
    tol::Float64 = 1e-6,
    panel_id::Union{Symbol, Nothing} = nothing,
    interior_beams::Vector{<:FrameElement} = FrameElement[]
)::Vector{TributaryLoad}
    isempty(shells) && return TributaryLoad[]
    isempty(beams) && return TributaryLoad[]
    
    # Extract panel boundary
    boundary = _get_shell_panel_boundary(shells; panel_id=panel_id)
    isempty(boundary) && error("Could not extract panel boundary from shells")
    
    # If interior beams specified, split polygon and compute sub-panel tributaries
    if !isempty(interior_beams)
        return _shell_to_tributary_loads_with_interior(
            boundary, beams, interior_beams, pressure;
            axis=axis, direction=direction, tol=tol
        )
    end
    
    # Convert boundary to Meshes.Point for tributary computation
    boundary_points = [Point(p[1]*u"m", p[2]*u"m") for p in boundary]
    
    # Compute tributary polygons
    trib_polys = isnothing(axis) ? 
        get_tributary_polygons(boundary_points) :
        get_tributary_polygons(boundary_points; axis=axis)
    
    # Map each beam to its edge
    beam_to_edge = Dict{FrameElement, Int}()
    unmatched_beams = FrameElement[]
    
    for beam in beams
        edge_idx = _match_beam_to_polygon_edge(beam, boundary_points; tol=tol)
        if isnothing(edge_idx)
            push!(unmatched_beams, beam)
        else
            beam_to_edge[beam] = edge_idx
        end
    end
    
    # Error if any beam doesn't match
    if !isempty(unmatched_beams)
        beam_ids = [b.id for b in unmatched_beams]
        error("Beams do not match any panel edge: $beam_ids. " *
              "Ensure all beams lie on the shell panel boundary, " *
              "or specify interior_beams for beams inside the panel.")
    end
    
    # Create TributaryLoad for each beam
    loads = TributaryLoad[]
    
    for (beam, edge_idx) in beam_to_edge
        tp = trib_polys[edge_idx]
        
        # Extract depth profile from polygon vertices
        positions, depths = _extract_depth_profile(tp)
        
        # Clamp negative depths to zero (can occur with numerical precision issues)
        depths = max.(depths, 0.0)
        widths = [w * u"m" for w in depths]
        
        trib_load = TributaryLoad(
            beam,
            positions,
            widths,
            pressure,
            direction
        )
        push!(loads, trib_load)
    end
    
    return loads
end

"""
Internal: Handle tributary computation with interior beams.
Splits the polygon along interior beam lines and computes tributary for each sub-panel.
Interior beams receive load from both sides.
"""
function _shell_to_tributary_loads_with_interior(
    boundary::Vector{NTuple{2, Float64}},
    beams::Vector{<:FrameElement},
    interior_beams::Vector{<:FrameElement},
    pressure::Quantity;
    axis::Union{Nothing, Vector{Float64}} = nothing,
    direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
    tol::Float64 = 1e-6
)::Vector{TributaryLoad}
    
    # Split the boundary polygon along interior beam lines
    sub_polygons = _split_polygon_by_beams(boundary, interior_beams; tol=tol)
    
    # For each sub-polygon, compute tributary to edge and interior beams
    # Interior beams will accumulate load from adjacent sub-polygons
    beam_loads = Dict{FrameElement, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}()
    
    for sub_poly in sub_polygons
        boundary_points = [Point(p[1]*u"m", p[2]*u"m") for p in sub_poly]
        
        # Compute tributary for this sub-panel
        trib_polys = isnothing(axis) ? 
            get_tributary_polygons(boundary_points) :
            get_tributary_polygons(boundary_points; axis=axis)
        
        # Match beams to sub-panel edges
        for beam in beams
            edge_idx = _match_beam_to_polygon_edge(beam, boundary_points; tol=tol)
            if !isnothing(edge_idx)
                tp = trib_polys[edge_idx]
                positions, depths = _extract_depth_profile(tp)
                depths = max.(depths, 0.0)
                
                if !haskey(beam_loads, beam)
                    beam_loads[beam] = Tuple{Vector{Float64}, Vector{Float64}}[]
                end
                push!(beam_loads[beam], (positions, depths))
            end
        end
    end
    
    # Create TributaryLoads, combining contributions for interior beams
    loads = TributaryLoad[]
    
    for (beam, contributions) in beam_loads
        if length(contributions) == 1
            # Edge beam - single contribution
            positions, depths = contributions[1]
            widths = [w * u"m" for w in depths]
            push!(loads, TributaryLoad(beam, positions, widths, pressure, direction))
        else
            # Interior beam - combine contributions from both sides
            combined = _combine_tributary_depths(contributions)
            widths = [w * u"m" for w in combined[2]]
            push!(loads, TributaryLoad(beam, combined[1], widths, pressure, direction))
        end
    end
    
    return loads
end
