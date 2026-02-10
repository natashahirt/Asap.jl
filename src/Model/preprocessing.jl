# =============================================================================
# DOF Index Population
# =============================================================================

"""
    populate_DOF_indices!(model::FrameModel)

Populate all indices, references, and other information in a frame model.
"""
function populate_DOF_indices!(model::FrameModel)
    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 6
    dofset = collect(0:n_dof - 1)

    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ dofset
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
    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 6
    dofset = collect(0:n_dof - 1)

    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ dofset
    end

    # Populate globalID for shell elements
    for element in model.elements
        populate_globalID!(element)
    end
end

"""
    populate_DOF_indices!(model::Model)

Populate all indices, references, and other information in a unified model.
"""
function populate_DOF_indices!(model::Model)
    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 6
    dofset = collect(0:n_dof - 1)

    # Assign global DOFs to nodes
    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ dofset
    end

    # Assign globalID to frame elements
    for element in model.frame_elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end

    # Assign globalID to shell elements
    for element in model.shell_elements
        populate_globalID!(element)
    end
end

"""
    populate_DOF_indices!(model::TrussModel)

Populate all indices, references, and other information in a truss model.
"""
function populate_DOF_indices!(model::TrussModel)
    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 3
    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof - 1)
    end

    for element in model.elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end
end

# =============================================================================
# Element Processing
# =============================================================================

"""
    process_elements!(model::FrameModel)
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a model.
"""
function process_elements!(model::FrameModel)
    for element in model.elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(model::ShellModel)
    
Process shell elements in a shell model.
"""
function process_elements!(model::ShellModel)
    for elem in model.elements
        process!(elem)
    end
end

"""
    process_elements!(model::Model)
    
Process all elements in a unified model (both frame and shell).
"""
function process_elements!(model::Model)
    # Process frame elements
    for element in model.frame_elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
    
    # Process shell elements
    for elem in model.shell_elements
        process!(elem)
    end
end

"""
    process_elements!(elements::Vector{<:FrameElement})
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a vector of elements.
"""
function process_elements!(elements::Vector{T}) where {T<:FrameElement}
    for element in elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(model::TrussModel)
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a truss model.
"""
function process_elements!(model::TrussModel)
    for element in model.elements
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(elements::Vector{<:ShellElement})

Process a vector of shell elements: compute LCS, transformation R, stiffness K, and mass M.
"""
function process_elements!(elements::Vector{T}) where {T<:ShellElement}
    for elem in elements
        process!(elem)
        populate_globalID!(elem)
    end
end

# =============================================================================
# Load Population
# =============================================================================

"""
    populate_load!(model::AbstractModel, load::NodeForce)

Populate the global load vector `model.P` with a nodal force.
"""
function populate_load!(model::AbstractModel, load::NodeForce)
    gid = load.node.globalID
    @inbounds for i in 1:3
        model.P[gid[i]] += to_newtons(load.value[i])
    end
end

"""
    populate_load!(P::Vector{Float64}, load::NodeForce)

Populate the global load vector P with a nodal force.
"""
function populate_load!(P::Vector{Float64}, load::NodeForce)
    gid = load.node.globalID
    @inbounds for i in 1:3
        P[gid[i]] += to_newtons(load.value[i])
    end
end

"""
    populate_load!(model, load::NodeMoment)

Populate the global load vector `model.P` with a nodal moment.
"""
function populate_load!(model::AbstractModel, load::NodeMoment)
    gid = load.node.globalID
    @inbounds for i in 1:3
        model.P[gid[i+3]] += to_newton_meters(load.value[i])
    end
end

"""
    populate_load!(P::Vector{Float64}, load::NodeMoment)

Populate the global load vector P with a nodal moment.
"""
function populate_load!(P::Vector{Float64}, load::NodeMoment)
    gid = load.node.globalID
    @inbounds for i in 1:3
        P[gid[i+3]] += to_newton_meters(load.value[i])
    end
end

"""
    populate_load!(model::FrameModel, load::ElementLoad)

Generate the fixed-end force vector `Q` for a given load, and populate the global fixed-end force vector `Pf`.
"""
function populate_load!(model::FrameModel, load::ElementLoad)
    idx = load.element.globalID
    fixed_end_forces_GCS = load.element.R' * q(load)
    load.element.Q += fixed_end_forces_GCS
    model.Pf[idx] += fixed_end_forces_GCS
end

"""
    populate_load!(model::Model, load::ElementLoad)

Generate the fixed-end force vector for a given load in a unified model.
"""
function populate_load!(model::Model, load::ElementLoad)
    idx = load.element.globalID
    fixed_end_forces_GCS = load.element.R' * q(load)
    load.element.Q += fixed_end_forces_GCS
    model.Pf[idx] += fixed_end_forces_GCS
end

"""
    populate_load!(Pf::Vector{Float64}, load::ElementLoad)

Populate an external fixed-end force vector Pf with respect to an elemental load.
"""
function populate_load!(Pf::Vector{Float64}, load::ElementLoad)
    idx = load.element.globalID
    Pf[idx] += load.element.R' * q(load)
end

"""
    populate_load!(model, load::AreaLoad)

Populate the global load vector with forces from an area load.

Behavior depends on `load.distribute_to`:
- `:nodes` → FEM approach: equivalent nodal forces on shell nodes
- `Vector{FrameElement}` → Tributary approach: compute tributary loads for beams
"""
function populate_load!(model::AbstractModel, load::AreaLoad)
    if load.distribute_to == :nodes
        # FEM approach: nodal forces on shell nodes
        for (node, fvec) in nodal_forces(load)
            idx = node.globalID[1:3]  # Translation DOFs
            model.P[idx] .+= fvec
        end
    else
        # Tributary approach: compute and apply tributary loads to beams
        if isnothing(load._tributary_loads)
            # Compute tributary loads (lazy evaluation)
            axis_vec = isnothing(load.axis) ? nothing : [load.axis[1], load.axis[2]]
            load._tributary_loads = _shell_to_tributary_loads(
                load.shells, 
                load.distribute_to, 
                load.pressure;
                axis=axis_vec,
                direction=load.direction,
                interior_beams=load.interior_beams
            )
        end
        
        # Apply each tributary load
        for trib_load in load._tributary_loads
            populate_load!(model, trib_load)
        end
    end
end

"""
    populate_load!(P::Vector{Float64}, load::AreaLoad)

Populate an external load vector with area load contributions.
"""
function populate_load!(P::Vector{Float64}, load::AreaLoad)
    if load.distribute_to == :nodes
        for (node, fvec) in nodal_forces(load)
            idx = node.globalID[1:3]
            P[idx] .+= fvec
        end
    else
        error("Tributary distribution requires model context. Use populate_load!(model, load) instead.")
    end
end

# =============================================================================
# Load Vector Assembly
# =============================================================================

"""
    populate_loads!(model::FrameModel)

Generate the nodal force vectors `model.P` (external) and `model.Pf` (fixed-end)
"""
function populate_loads!(model::FrameModel)
    fill!(model.P, 0.0)
    fill!(model.Pf, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::ShellModel)

Generate the nodal force vectors `model.P` for a shell model.
"""
function populate_loads!(model::ShellModel)
    fill!(model.P, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::Model)

Generate the nodal force vectors for a unified model.
"""
function populate_loads!(model::Model)
    fill!(model.P, 0.0)
    fill!(model.Pf, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::TrussModel)

Generate the nodal force vectors `model.P`
"""
function populate_loads!(model::TrussModel)
    fill!(model.P, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

# =============================================================================
# Load Vector Creation (for external load vectors)
# =============================================================================

"""
    create_F(model::FrameModel, loads::Vector{AbstractLoad})

Create load vector F = P - Pf for a given load set.
"""
function create_F(model::FrameModel, loads::Vector{T}) where {T<:AbstractLoad}
    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if load isa ElementLoad
            populate_load!(Pf, load)
        else
            populate_load!(P, load)
        end
    end

    return P - Pf
end

"""
    create_F(model::Model, loads::Vector{AbstractLoad})

Create load vector F = P - Pf for a given load set.
"""
function create_F(model::Model, loads::Vector{T}) where {T<:AbstractLoad}
    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if load isa ElementLoad
            populate_load!(Pf, load)
        else
            populate_load!(P, load)
        end
    end

    return P - Pf
end

"""
    create_F(model::TrussModel, loads::Vector{NodeForce})

Create load vector F = P for truss models.
"""
function create_F(model::TrussModel, loads::Vector{NodeForce})
    P = zeros(model.nDOFs)

    for load in loads
        populate_load!(P, load)
    end

    return P
end

# =============================================================================
# Stiffness Matrix Assembly
# =============================================================================

"""
    create_S!(model::FrameModel)

Assemble the global stiffness matrix S for a frame model.
"""
function create_S!(model::FrameModel)
    nnz_tot = 144 * length(model.elements)  # 12×12 per frame element (exact)
    I = Vector{Int64}(undef, nnz_tot)
    J = Vector{Int64}(undef, nnz_tot)
    V = Vector{Float64}(undef, nnz_tot)

    pos = 0
    for element in model.elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i = 1:n, j = 1:n
            pos += 1
            I[pos] = idx[i]
            J[pos] = idx[j]
            V[pos] = element.K[i,j]
        end
    end

    model.S = sparse(I, J, V, model.nDOFs, model.nDOFs)
end

"""
    create_S!(model::ShellModel)

Assemble the global stiffness matrix S for a shell model.
"""
function create_S!(model::ShellModel)
    nnz_tot = sum(length(e.globalID)^2 for e in model.elements; init=0)
    I = Vector{Int64}(undef, nnz_tot)
    J = Vector{Int64}(undef, nnz_tot)
    V = Vector{Float64}(undef, nnz_tot)

    pos = 0
    for element in model.elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i = 1:n, j = 1:n
            pos += 1
            I[pos] = idx[i]
            J[pos] = idx[j]
            V[pos] = element.K[i,j]
        end
    end

    model.S = sparse(I, J, V, model.nDOFs, model.nDOFs)
end

"""
    create_S!(model::Model)

Assemble the global stiffness matrix S for a unified model.
"""
function create_S!(model::Model)
    nnz_tot = 144 * length(model.frame_elements) +
              sum(length(e.globalID)^2 for e in model.shell_elements; init=0)
    I = Vector{Int64}(undef, nnz_tot)
    J = Vector{Int64}(undef, nnz_tot)
    V = Vector{Float64}(undef, nnz_tot)

    pos = 0
    for element in model.frame_elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i = 1:n, j = 1:n
            pos += 1
            I[pos] = idx[i]
            J[pos] = idx[j]
            V[pos] = element.K[i,j]
        end
    end

    for element in model.shell_elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i = 1:n, j = 1:n
            pos += 1
            I[pos] = idx[i]
            J[pos] = idx[j]
            V[pos] = element.K[i,j]
        end
    end

    model.S = sparse(I, J, V, model.nDOFs, model.nDOFs)
end

"""
    create_S!(model::TrussModel)

Assemble the global stiffness matrix S for a truss model.
"""
function create_S!(model::TrussModel)
    nnz_tot = 36 * length(model.elements)  # 6×6 per truss element (exact)
    I = Vector{Int64}(undef, nnz_tot)
    J = Vector{Int64}(undef, nnz_tot)
    V = Vector{Float64}(undef, nnz_tot)

    pos = 0
    for element in model.elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i = 1:n, j = 1:n
            pos += 1
            I[pos] = idx[i]
            J[pos] = idx[j]
            V[pos] = element.K[i,j]
        end
    end

    model.S = sparse(I, J, V, model.nDOFs, model.nDOFs)
end

# =============================================================================
# General Assembly Functions (for direct element vectors)
# =============================================================================

"""
    assemble_stiffness(elements, n_dof) -> SparseMatrixCSC

Assemble global stiffness matrix from a vector of elements.
Works with any element type that has `globalID` and `K` fields.
"""
function assemble_stiffness(elements::Vector, n_dof::Int)
    nnz_est = sum(length(e.globalID)^2 for e in elements; init=0)
    I = Int[]; sizehint!(I, nnz_est)
    J = Int[]; sizehint!(J, nnz_est)
    V = Float64[]; sizehint!(V, nnz_est)
    
    for elem in elements
        idx = elem.globalID
        n = length(idx)
        
        for i = 1:n, j = 1:n
            push!(I, idx[i])
            push!(J, idx[j])
            push!(V, elem.K[i, j])
        end
    end
    
    return sparse(I, J, V, n_dof, n_dof)
end

"""
    assemble_stiffness!(S, elements)

Add element stiffness contributions to existing sparse matrix S.
"""
function assemble_stiffness!(S::SparseMatrixCSC, elements::Vector)
    for elem in elements
        idx = elem.globalID
        n = length(idx)
        
        for i = 1:n, j = 1:n
            S[idx[i], idx[j]] += elem.K[i, j]
        end
    end
    
    return S
end

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
function add_springs!(model::AbstractModel, springs::Vector{Spring})
    # Populate spring IDs and global DOF indices
    make_ids!(springs)
    for spring in springs
        populate_globalID!(spring)
    end
    
    # Add spring stiffness to global matrix
    assemble_stiffness!(model.S, springs)
    
    return model
end

"""
    add_springs!(model::AbstractModel, spring::Spring)

Add a single spring to a model. Convenience wrapper.
"""
add_springs!(model::AbstractModel, spring::Spring) = add_springs!(model, [spring])