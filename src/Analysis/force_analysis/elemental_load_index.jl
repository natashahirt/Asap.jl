# =============================================================================
# ElementInternalForces — sizehinted batch constructors
# =============================================================================


"""
    get_elemental_loads(model) -> Vector{Vector{AbstractLoad}}

Return a per-element load map.  The result is lazily cached on the model
(`model._elemental_loads`) and invalidated automatically by `process!`.
"""
function get_elemental_loads(model::FrameModel)
    cache = model._elemental_loads
    cache !== nothing && return cache
    element_to_loads = [AbstractLoad[] for _ in 1:model.nElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        push!(element_to_loads[load.element.elementID], load)
    end
    model._elemental_loads = element_to_loads
    return element_to_loads
end

function get_elemental_loads(model::Model)
    cache = model._elemental_loads
    cache !== nothing && return cache
    element_to_loads = [AbstractLoad[] for _ in 1:model.nFrameElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        if load.element isa FrameElement
            push!(element_to_loads[load.element.elementID], load)
        end
    end
    model._elemental_loads = element_to_loads
    return element_to_loads
end

"""
    get_elemental_loads(loads, n_elements) -> Vector{Vector{AbstractLoad}}

Build a per-element load map from an arbitrary load vector (no model required).
Useful for post-processing displacement vectors from `solve(model, cases)` where
each case has its own load set.
"""
function get_elemental_loads(loads::Vector{<:AbstractLoad}, n_elements::Int)
    element_to_loads = [AbstractLoad[] for _ in 1:n_elements]
    for load in loads
        hasproperty(load, :element) || continue
        eid = load.element.elementID
        (1 <= eid <= n_elements) && push!(element_to_loads[eid], load)
    end
    return element_to_loads
end

"""Invalidate the cached elemental loads map (call when loads or elements change)."""
function clear_elemental_loads!(model::AbstractModel)
    hasproperty(model, :_elemental_loads) && (model._elemental_loads = nothing)
end

