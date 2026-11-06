# Per-thread buffer for populate_load! (avoids R' * q allocation per load).
# Filled in `Asap.__init__` (see threading.jl).
const _POPULATE_LOAD_BUF = Vector{Float64}[]

# =============================================================================
# DOF Index Population
# =============================================================================

"""Populate model.DOFs as a flat Bool vector from node.dof — zero extra allocation."""
function _populate_DOFs_flat!(model, n_nodes::Int, n_dof_per_node::Int)
    total = n_nodes * n_dof_per_node
    DOFs = Vector{Bool}(undef, total)
    @inbounds for (i, node) in enumerate(model.nodes)
        base = (i - 1) * n_dof_per_node
        dofs = node.dof
        for j in 1:n_dof_per_node
            DOFs[base + j] = dofs[j]
        end
    end
    model.DOFs = DOFs
    model.nDOFs = total
    model.nNodes = n_nodes
    model.freeDOFs = findall(DOFs)
    model.fixedDOFs = findall(.!DOFs)

    # `process!` is also the topology-change boundary. Keep all nodal state
    # vectors synchronized when normalization adds or removes nodes.
    length(model.P) == total || (model.P = zeros(total))
    hasproperty(model, :Pf) &&
        length(model.Pf) != total && (model.Pf = zeros(total))
    length(model.u) == total || (model.u = zeros(total))
    length(model.reactions) == total || (model.reactions = zeros(total))
end

"""
    populate_DOF_indices!(model::FrameModel)

Populate all indices, references, and other information in a frame model.
"""
function populate_DOF_indices!(model::FrameModel)
    n_dof = 6
    nn = length(model.nodes)
    _populate_DOFs_flat!(model, nn, n_dof)
    model.nElements = length(model.elements)

    for (i, node) in enumerate(model.nodes)
        base = (i - 1) * n_dof
        node.globalID = [base + 1, base + 2, base + 3, base + 4, base + 5, base + 6]
    end

    for element in model.elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end
end

"""
    populate_DOF_indices!(model::ShellModel)

Populate all indices, references, and other information in a shell model.
"""
function populate_DOF_indices!(model::ShellModel)
    n_dof = 6
    nn = length(model.nodes)
    _populate_DOFs_flat!(model, nn, n_dof)
    model.nElements = length(model.elements)

    for (i, node) in enumerate(model.nodes)
        base = (i - 1) * n_dof
        node.globalID = [base + 1, base + 2, base + 3, base + 4, base + 5, base + 6]
    end

    for element in model.elements
        populate_globalID!(element)
    end
end

"""
    populate_DOF_indices!(model::Model)

Populate all indices, references, and other information in a unified model.
"""
function populate_DOF_indices!(model::Model)
    n_dof = 6
    nn = length(model.nodes)
    _populate_DOFs_flat!(model, nn, n_dof)

    for (i, node) in enumerate(model.nodes)
        base = (i - 1) * n_dof
        node.globalID = [base + 1, base + 2, base + 3, base + 4, base + 5, base + 6]
    end

    for element in model.frame_elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end

    for element in model.shell_elements
        populate_globalID!(element)
    end
end

"""
    populate_DOF_indices!(model::TrussModel)

Populate all indices, references, and other information in a truss model.
"""
function populate_DOF_indices!(model::TrussModel)
    n_dof = 3
    nn = length(model.nodes)
    _populate_DOFs_flat!(model, nn, n_dof)
    model.nElements = length(model.elements)

    for (i, node) in enumerate(model.nodes)
        base = (i - 1) * n_dof
        node.globalID = [base + 1, base + 2, base + 3]
    end

    for element in model.elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end
end

