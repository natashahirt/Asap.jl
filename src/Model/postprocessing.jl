using LinearAlgebra: mul!

# =============================================================================
# Reactions
# =============================================================================

"""
    reactions!(model::AbstractModel)

Populate external reaction forces in `model.reactions`
"""
function reactions!(model::FrameModel)
    if length(model.reactions) == model.nDOFs
        fill!(model.reactions, 0.0)
    else
        model.reactions = zeros(model.nDOFs)
    end
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u + model.Pf[model.fixedDOFs]
end

function reactions!(model::ShellModel)
    if length(model.reactions) == model.nDOFs
        fill!(model.reactions, 0.0)
    else
        model.reactions = zeros(model.nDOFs)
    end
    # Subtract applied forces at supports: true reaction = K·u − P at fixed DOFs
    # (AreaLoad distributes pressure to ALL nodes, including supports)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u - model.P[model.fixedDOFs]
end

function reactions!(model::Model)
    if length(model.reactions) == model.nDOFs
        fill!(model.reactions, 0.0)
    else
        model.reactions = zeros(model.nDOFs)
    end
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u + model.Pf[model.fixedDOFs]
end

function reactions!(model::TrussModel)
    if length(model.reactions) == model.nDOFs
        fill!(model.reactions, 0.0)
    else
        model.reactions = zeros(model.nDOFs)
    end
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u
end

# =============================================================================
# Node Post-Processing
# =============================================================================

"""
    _post_process_nodes_6dof!(model)

Shared 6-DOF node post-processing (Frame, Shell, Model).
Populate nodal reaction and displacement fields with Unitful quantities.
Reactions are in SI base units: forces in N, moments in N*m.
Displacements are in SI base units: translations in m, rotations in radians.
"""
function _post_process_nodes_6dof!(model)
    rxn = model.reactions
    u_vec = model.u
    for node in model.nodes
        gid = node.globalID
        n = length(gid)
        rq = Vector{Quantity}(undef, n)
        dq = Vector{Quantity}(undef, n)
        @inbounds for i in 1:n
            g = gid[i]
            if i <= 3
                rq[i] = rxn[g] * u"N"
                dq[i] = u_vec[g] * u"m"
            else
                rq[i] = rxn[g] * u"N*m"
                dq[i] = u_vec[g] * u"rad"
            end
        end
        node.reaction = rq
        node.displacement = dq
    end
end

post_process_nodes!(model::FrameModel) = _post_process_nodes_6dof!(model)
post_process_nodes!(model::ShellModel) = _post_process_nodes_6dof!(model)
post_process_nodes!(model::Model) = _post_process_nodes_6dof!(model)

function post_process_nodes!(model::TrussModel)
    rxn = model.reactions
    u_vec = model.u
    for node in model.nodes
        gid = node.globalID
        n = length(gid)
        rq = Vector{Quantity}(undef, n)
        dq = Vector{Quantity}(undef, n)
        @inbounds for i in 1:n
            g = gid[i]
            rq[i] = rxn[g] * u"N"
            dq[i] = u_vec[g] * u"m"
        end
        node.reaction = rq
        node.displacement = dq
    end
end

# =============================================================================
# Element Post-Processing
# =============================================================================

"""
    post_process_elements!(model::AbstractModel)

Populate elemental LCS force vectors in `element.forces`.
Uses preallocated buffers to avoid per-element allocations.
"""
function post_process_elements!(model::FrameModel)
    u_e  = Vector{Float64}(undef, 12)  # Element displacement buffer
    Ku_e = Vector{Float64}(undef, 12)  # K * u buffer
    @inbounds for element in model.elements
        idx = element.globalID
        for i in 1:12; u_e[i] = model.u[idx[i]]; end
        mul!(Ku_e, element.K, u_e)
        for i in 1:12; Ku_e[i] += element.Q[i]; end
        mul!(element.forces, element.R, Ku_e)
    end
end

function post_process_elements!(model::ShellModel)
    # Shell elements don't have .forces field in same way
    # Stresses are computed via stress() function instead
end

function post_process_elements!(model::Model)
    # Frame elements
    u_e  = Vector{Float64}(undef, 12)
    Ku_e = Vector{Float64}(undef, 12)
    @inbounds for element in model.frame_elements
        idx = element.globalID
        for i in 1:12; u_e[i] = model.u[idx[i]]; end
        mul!(Ku_e, element.K, u_e)
        for i in 1:12; Ku_e[i] += element.Q[i]; end
        mul!(element.forces, element.R, Ku_e)
    end
    
    # Shell elements - stresses computed via stress() function
end

function post_process_elements!(model::TrussModel)
    n_dof = 6  # 3 DOFs per node × 2 nodes for truss
    u_e  = Vector{Float64}(undef, n_dof)
    Ku_e = Vector{Float64}(undef, n_dof)
    @inbounds for element in model.elements
        idx = element.globalID
        for i in 1:n_dof; u_e[i] = model.u[idx[i]]; end
        mul!(Ku_e, element.K, u_e)
        mul!(element.forces, element.R, Ku_e)
    end
end

# =============================================================================
# Main Post-Processing Entry Point
# =============================================================================

"""
    post_process!(model::AbstractModel; targets::Symbol = :all)

Post process a model after solving for displacements using `solve!(model)`.

# Target Symbols
- `:all`      — reactions + node Unitful fields + element forces (default)
- `:elements` — element forces only (fastest for internal analysis like EFM)
- `:nodes`    — reactions + node Unitful fields only
- `:none`     — skip all post-processing
"""
function post_process!(model::AbstractModel; targets::Symbol = :all)
    if targets ∈ (:all, :nodes)
        reactions!(model)
        post_process_nodes!(model)
    end
    if targets ∈ (:all, :elements)
        post_process_elements!(model)
    end
end
