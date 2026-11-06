# =============================================================================
# Spring Element Support
# =============================================================================

"""
    populate_globalID!(spring::Spring)

Set the global DOF indices for a spring element based on its attached node.
"""
function populate_globalID!(spring::Spring)
    spring.globalID = spring.node.globalID
end

"""
    make_ids!(springs::Vector{Spring})

Assign element IDs to springs.
"""
function make_ids!(springs::Vector{Spring})
    for (i, spring) in enumerate(springs)
        spring.elementID = i
    end
end

"""
    add_springs!(model::AbstractModel, springs::Vector{Spring})

Add grounded springs to a model after processing.

Call this after `process!(model)` but before `solve!(model)`.
Springs add diagonal stiffness terms at their attached nodes' DOFs.

# Example
```julia
model = Model(nodes, elements, loads)
process!(model)

# Add springs
springs = [Spring(node1, 1000u"kN/m"), Spring(node2; kz=500u"kN/m")]
add_springs!(model, springs)

solve!(model)
```
"""
function _assemble_model_springs!(model::Union{FrameModel,Model})
    isempty(model.springs) && return model
    make_ids!(model.springs)
    for spring in model.springs
        populate_globalID!(spring)
    end
    assemble_stiffness!(model.S, model.springs)
    model._scatter = nothing
    clear_factorization!(model)
    return model
end

function add_springs!(
    model::Union{FrameModel,Model},
    springs::Vector{Spring},
)
    # Remember each spring once. Stiffness is assembled on every call: a
    # `values_only` update rebuilds `S` from elements and drops the spring
    # diagonals, and the no-tension Winkler loop relies on this call to put
    # the active subset back. Skipping assembly when the object is already
    # stored leaves that solve with no soil stiffness.
    for spring in springs
        if !any(existing -> existing === spring, model.springs)
            push!(model.springs, spring)
        end
    end
    model.processed || return model

    make_ids!(springs)
    for spring in springs
        populate_globalID!(spring)
    end
    assemble_stiffness!(model.S, springs)
    model._scatter = nothing
    clear_factorization!(model)
    return model
end

function add_springs!(
    model::Union{ShellModel,TrussModel},
    springs::Vector{Spring},
)
    make_ids!(springs)
    for spring in springs
        populate_globalID!(spring)
    end
    assemble_stiffness!(model.S, springs)
    model._scatter = nothing
    clear_factorization!(model)
    return model
end

"""
    add_springs!(model::AbstractModel, spring::Spring)

Add a single spring to a model. Convenience wrapper.
"""
add_springs!(
    model::AbstractModel,
    spring::Spring,
) = add_springs!(model, [spring])