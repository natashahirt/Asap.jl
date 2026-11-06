# =============================================================================
# _update_S_values!  —  in-place stiffness scatter (same sparsity pattern)
# =============================================================================

"""Recompute every shell `global_K!` on a concretely typed element vector."""
function _recompute_shell_K!(shells::Vector{Sh}) where {Sh}
    Threads.@threads for i in eachindex(shells)
        global_K!(shells[i])
    end
    return shells
end

"""
    _update_S_values!(model)

Zero the existing sparse matrix `model.S` and scatter updated element K
values back into it.  Reuses the CSC sparsity pattern from the initial
`create_S!` call — avoids COO allocation + `sparse()` reconstruction.

Uses the cached [`_ScatterMap`](@ref) from `create_S!` when it still matches
`model.S`; otherwise falls back to indexed `S[i, j] +=` writes. Both paths add
the element contributions in the same order.
"""
# The element vectors are stored in abstractly typed fields
# (`Vector{<:FrameElement}`, `Vector{<:ShellElement}`). A barrier on the
# runtime vector type lets the scatter loop see a concrete `K`.
_update_S_values!(model::Union{FrameModel, ShellModel, TrussModel}) =
    _assemble_S!(model, model.elements)

_update_S_values!(model::Model) =
    _assemble_S!(model, model.frame_elements, model.shell_elements)

function _assemble_S!(model::AbstractModel, elements::Vector{E}) where {E}
    S = model.S
    fill!(S.nzval, 0.0)
    sm = model._scatter
    n_expected = _coo_length(elements)
    if sm !== nothing && sm.S === S && sm.nnz == nnz(S) && length(sm.idx) == n_expected
        _scatter_values!(S.nzval, sm.idx, 0, elements)
    else
        _scatter_values_indexed!(S, elements)
    end
    return S
end

function _assemble_S!(
    model::AbstractModel, frames::Vector{F}, shells::Vector{Sh},
) where {F, Sh}
    S = model.S
    fill!(S.nzval, 0.0)
    sm = model._scatter
    n_expected = _coo_length(frames) + _coo_length(shells)
    if sm !== nothing && sm.S === S && sm.nnz == nnz(S) && length(sm.idx) == n_expected
        # Frames first, then shells: the same order `create_S!` used for the
        # mixed model, so the scatter-map slots stay aligned.
        p = _scatter_values!(S.nzval, sm.idx, 0, frames)
        _scatter_values!(S.nzval, sm.idx, p, shells)
    else
        _scatter_values_indexed!(S, frames)
        _scatter_values_indexed!(S, shells)
    end
    return S
end

_coo_length(elements::Vector) = sum(e -> length(e.globalID)^2, elements; init = 0)

function _scatter_values!(
    nzval::Vector{Float64}, idx::Vector{Int}, p::Int, elements::Vector{E},
) where {E}
    @inbounds for element in elements
        K = element.K
        n = length(element.globalID)
        for i in 1:n, j in 1:n
            p += 1
            nzval[idx[p]] += K[i, j]
        end
    end
    return p
end

function _scatter_values_indexed!(S::SparseMatrixCSC, elements::Vector{E}) where {E}
    for element in elements
        idx = element.globalID
        K = element.K
        n = length(idx)
        @inbounds for i in 1:n, j in 1:n
            S[idx[i], idx[j]] += K[i, j]
        end
    end
end
