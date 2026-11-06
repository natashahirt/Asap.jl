
# =============================================================================
# Shell Panel Boundary Extraction (Internal)
# =============================================================================

"""Internal: Extract outer boundary polygon from shell elements."""
function _get_shell_panel_boundary(
    shells::Vector{<:ShellElement};
    panel_id::Union{Symbol, Nothing} = nothing
)::Vector{NTuple{2, Float64}}
    # Filter by panel_id if specified
    filtered = isnothing(panel_id) ? shells : filter(s -> s.id == panel_id, shells)
    isempty(filtered) && return NTuple{2, Float64}[]
    
    # Collect all edges (as sorted node position pairs to identify shared edges)
    # Edge format: ((x1,y1), (x2,y2)) with (x1,y1) < (x2,y2) lexicographically
    edge_count = Dict{Tuple{NTuple{2,Float64}, NTuple{2,Float64}}, Int}()
    edge_to_original = Dict{Tuple{NTuple{2,Float64}, NTuple{2,Float64}}, NTuple{2, NTuple{2,Float64}}}()
    
    for shell in filtered
        nodes = shell.nodes
        n_nodes = length(nodes)
        
        for i in 1:n_nodes
            j = mod1(i + 1, n_nodes)
            n1, n2 = nodes[i], nodes[j]
            
            # Node stores position as a vector: position[1]=x, position[2]=y, position[3]=z
            p1 = (Float64(ustrip(u"m", n1.position[1])), Float64(ustrip(u"m", n1.position[2])))
            p2 = (Float64(ustrip(u"m", n2.position[1])), Float64(ustrip(u"m", n2.position[2])))
            
            # Create canonical edge key (sorted)
            key = p1 < p2 ? (p1, p2) : (p2, p1)
            original = (p1, p2)  # Keep original orientation
            
            edge_count[key] = get(edge_count, key, 0) + 1
            if !haskey(edge_to_original, key)
                edge_to_original[key] = original
            end
        end
    end
    
    # Keep only boundary edges (count == 1)
    boundary_edges = [(edge_to_original[k], k) for (k, count) in edge_count if count == 1]
    isempty(boundary_edges) && return NTuple{2, Float64}[]
    
    # Chain edges into a polygon
    # Build adjacency: point → list of edges starting/ending there
    point_to_edges = Dict{NTuple{2,Float64}, Vector{Int}}()
    edges = [e[1] for e in boundary_edges]  # Original orientation
    
    for (idx, (p1, p2)) in enumerate(edges)
        push!(get!(point_to_edges, p1, Int[]), idx)
        push!(get!(point_to_edges, p2, Int[]), idx)
    end
    
    # Walk the boundary starting from first edge
    polygon = NTuple{2, Float64}[]
    used = Set{Int}()
    
    current_edge_idx = 1
    current_point = edges[1][1]
    push!(polygon, current_point)
    push!(used, 1)
    current_point = edges[1][2]
    push!(polygon, current_point)
    
    while length(used) < length(edges)
        # Find next edge that shares current_point
        found = false
        for edge_idx in get(point_to_edges, current_point, Int[])
            if edge_idx ∉ used
                push!(used, edge_idx)
                p1, p2 = edges[edge_idx]
                # Move to the other end of this edge
                next_point = (p1 == current_point) ? p2 : p1
                push!(polygon, next_point)
                current_point = next_point
                found = true
                break
            end
        end
        !found && break
    end
    
    # Remove duplicate last point if it matches first
    if length(polygon) > 1 && polygon[end] == polygon[1]
        pop!(polygon)
    end
    
    # Ensure CCW orientation
    if _polygon_area(polygon) < 0
        reverse!(polygon)
    end
    
    return polygon
end
