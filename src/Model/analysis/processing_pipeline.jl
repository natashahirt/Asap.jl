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
    model.processed = false
    normalize_parametric_topology!(model)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    make_ids!(model)
    process_elements!(model)

    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)
    _assemble_model_springs!(model)

    model._factorization = nothing
    model._elemental_loads = nothing
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

    model._factorization = nothing
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
    model.processed = false
    normalize_parametric_topology!(model)
    model.nNodes = length(model.nodes)
    model.nFrameElements = length(model.frame_elements)
    model.nShellElements = length(model.shell_elements)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)
    _assemble_model_springs!(model)

    model._factorization = nothing
    model._elemental_loads = nothing
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

    model._factorization = nothing
    model.processed = true
end

