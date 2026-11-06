# =============================================================================
# Stiffness Matrix Assembly
# =============================================================================

"""
    _set_S!(model, I, J, V)

Assemble `model.S` from COO triplets and cache the `nzval` position of every
triplet so later value-only updates can skip the sparse search.
"""
function _set_S!(model::AbstractModel, I::Vector{Int64}, J::Vector{Int64}, V::Vector{Float64})
    S = sparse(I, J, V, model.nDOFs, model.nDOFs)
    model.S = S
    model._scatter = _build_scatter_map(S, I, J)
    return S
end

"""
    _build_scatter_map(S, I, J) -> Union{_ScatterMap, Nothing}

Locate each triplet `(I[p], J[p])` in the CSC structure of `S`. Returns `nothing`
if any triplet has no stored slot (cannot happen for `sparse(I, J, V)`, which
keeps numerical zeros, but the in-place path must never write to a wrong slot).
"""
function _build_scatter_map(S::SparseMatrixCSC{Float64,Int64}, I::Vector{Int64}, J::Vector{Int64})
    idx = Vector{Int}(undef, length(I))
    colptr, rowval = S.colptr, S.rowval
    ok = Threads.Atomic{Bool}(true)
    Threads.@threads for p in eachindex(I)
        @inbounds begin
            lo, hi = colptr[J[p]], colptr[J[p] + 1] - 1
            k = searchsortedfirst(view(rowval, lo:hi), I[p]) + lo - 1
            if k <= hi && rowval[k] == I[p]
                idx[p] = k
            else
                ok[] = false
            end
        end
    end
    return ok[] ? _ScatterMap(S, nnz(S), idx) : nothing
end

