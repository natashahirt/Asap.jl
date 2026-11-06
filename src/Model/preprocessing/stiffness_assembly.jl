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

    _set_S!(model, I, J, V)
end

"""
    create_S!(model::ShellModel)

Assemble the global stiffness matrix S for a shell model.
Thread-parallel: each element fills its own contiguous block of COO triplets.
"""
function create_S!(model::ShellModel)
    n_elem = length(model.elements)

    # Precompute per-element nnz and cumulative offsets
    elem_nnz = Vector{Int}(undef, n_elem)
    @inbounds for k in 1:n_elem
        elem_nnz[k] = length(model.elements[k].globalID)^2
    end
    offsets = Vector{Int}(undef, n_elem + 1)
    offsets[1] = 0
    @inbounds for k in 1:n_elem
        offsets[k + 1] = offsets[k] + elem_nnz[k]
    end
    nnz_tot = offsets[end]

    I = Vector{Int64}(undef, nnz_tot)
    J = Vector{Int64}(undef, nnz_tot)
    V = Vector{Float64}(undef, nnz_tot)

    Threads.@threads for k in 1:n_elem
        element = model.elements[k]
        idx = element.globalID
        n = length(idx)
        p = offsets[k]
        @inbounds for i = 1:n, j = 1:n
            p += 1
            I[p] = idx[i]
            J[p] = idx[j]
            V[p] = element.K[i,j]
        end
    end

    _set_S!(model, I, J, V)
end

"""
    create_S!(model::Model)

Assemble the global stiffness matrix S for a unified model.
Frame elements are filled sequentially (few elements), then shell elements
are filled in parallel (each element writes to its own contiguous COO block).
"""
function create_S!(model::Model)
    n_frame = length(model.frame_elements)
    n_shell = length(model.shell_elements)
    nnz_frame = 144 * n_frame

    # Precompute per-shell-element nnz and cumulative offsets
    shell_nnz = Vector{Int}(undef, n_shell)
    @inbounds for k in 1:n_shell
        shell_nnz[k] = length(model.shell_elements[k].globalID)^2
    end
    shell_offsets = Vector{Int}(undef, n_shell + 1)
    shell_offsets[1] = nnz_frame
    @inbounds for k in 1:n_shell
        shell_offsets[k + 1] = shell_offsets[k] + shell_nnz[k]
    end
    nnz_tot = shell_offsets[end]

    I = Vector{Int64}(undef, nnz_tot)
    J = Vector{Int64}(undef, nnz_tot)
    V = Vector{Float64}(undef, nnz_tot)

    # Frame elements — sequential (typically few)
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

    # Shell elements — parallel COO fill
    Threads.@threads for k in 1:n_shell
        element = model.shell_elements[k]
        idx = element.globalID
        n = length(idx)
        p = shell_offsets[k]
        @inbounds for i = 1:n, j = 1:n
            p += 1
            I[p] = idx[i]
            J[p] = idx[j]
            V[p] = element.K[i,j]
        end
    end

    _set_S!(model, I, J, V)
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

    _set_S!(model, I, J, V)
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

