# =============================================================================
# Model Processing
# =============================================================================

"""
    process!(model::FrameModel)

Process a structural model: add linkages between nodes and elements, 
determine DOF orders, generate the load vectors P, Pf, and assemble 
the global stiffness matrix, S.
"""
function process!(model::FrameModel)
    make_ids!(model)

    if any(typeof.(model.elements) .<: BridgeElement)
        processBridge!(model)
        make_ids!(model)
    else
        process_elements!(model)
    end

    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model.processed = true
end

"""
    process!(model::ShellModel)

Process a shell model: add linkages between nodes and elements,
determine DOF orders, generate the load vector P, and assemble
the global stiffness matrix, S.
"""
function process!(model::ShellModel)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model.processed = true
end

"""
    process!(model::Model)

Process a unified model: add linkages between nodes and elements,
determine DOF orders, generate the load vectors P, Pf, and assemble
the global stiffness matrix, S.

Works with mixed frame+shell models or single-type models.
"""
function process!(model::Model)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model.processed = true
end

"""
    process!(model::TrussModel)

Process a structural truss model.
"""
function process!(model::TrussModel)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model.processed = true
end

# =============================================================================
# Solving
# =============================================================================

"""
    solve!(model::FrameModel; reprocess = false)

Solve for the nodal displacements of a structural model. 
`reprocess = true` reevaluates all node/element properties and reassembles the global stiffness matrix.
"""
function solve!(model::FrameModel; reprocess = false)
    if !model.processed || reprocess
        for element in model.elements
            if typeof(element) <: Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    model.compliance = U' * F
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

"""
    solve!(model::ShellModel; reprocess = false)

Solve for the nodal displacements of a shell model.
"""
function solve!(model::ShellModel; reprocess = false)
    if !model.processed || reprocess
        process!(model)
    end

    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    model.compliance = U' * model.P[idx]
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

"""
    solve!(model::Model; reprocess = false)

Solve for the nodal displacements of a unified model.
Works with mixed frame+shell models or single-type models.
"""
function solve!(model::Model; reprocess = false)
    if !model.processed || reprocess
        # Clear fixed-end forces for frame elements
        for element in model.frame_elements
            if typeof(element) <: Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    model.compliance = U' * F
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

"""
    solve!(model::TrussModel; reprocess = false)

Solve for the nodal displacements of a structural truss model.
"""
function solve!(model::TrussModel; reprocess = false)
    if !model.processed || reprocess
        process!(model)
    end

    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    model.compliance = U' * model.P[idx]
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

# =============================================================================
# Solving with Custom Loads
# =============================================================================

"""
    solve(model::FrameModel, L::Vector{AbstractLoad})

Return the displacement vector under a given load set L.
"""
function solve(model::FrameModel, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve(model::Model, L::Vector{AbstractLoad})

Return the displacement vector under a given load set L.
"""
function solve(model::Model, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::FrameModel, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::FrameModel, L::Vector{<:AbstractLoad})
    model.loads = L
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve!(model::Model, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::Model, L::Vector{<:AbstractLoad})
    model.loads = L
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector under a given load set L.
"""
function solve(model::TrussModel, L::Vector{NodeForce})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::TrussModel, L::Vector{NodeForce})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::TrussModel, L::Vector{NodeForce})
    model.loads = L
    process!(model)
    solve!(model)
    post_process!(model)
end
