# =============================================================================
# Factorization Cache Helpers
# =============================================================================

"""
    factorize!(model::AbstractModel)

Compute and cache the stiffness factorization for the free DOFs.
Reuse with `solve_factorized!` to avoid re-factorizing when only loads change.

Uses Cholesky (CHOLMOD) for SPD matrices, LDLT (CHOLMOD) for symmetric
indefinite, falls back to LU (UMFPACK).
"""
function factorize!(model::AbstractModel)
    model.processed || error("Model must be processed before factorizing. Call process!(model) first.")
    idx = model.freeDOFs
    log_analysis!(model, :factor; n_free = length(idx))
    K = Symmetric(model.S[idx, idx])
    fact = cholesky(K; check=false)
    if !issuccess(fact)
        @warn "Stiffness matrix not SPD — trying LDLᵀ"
        fact = ldlt(K; check=false)
        if !issuccess(fact)
            @warn "LDLᵀ failed — falling back to LU"
            fact = lu(model.S[idx, idx])
        end
    end
    model._factorization = fact
    return model._factorization
end

"""
    clear_factorization!(model::AbstractModel)

Invalidate the cached factorization (call after stiffness changes).
"""
function clear_factorization!(model::AbstractModel)
    model._factorization = nothing
    model._factorization_stale = false
    model._factor_key = UInt(0)
    return model
end

# =============================================================================
# Internal Analysis Implementations
# =============================================================================

# Helper: get or compute factorization for free DOFs.
# Cholesky (CHOLMOD) is fastest for symmetric positive-definite stiffness
# matrices — the common case.  Falls back to LDLᵀ (CHOLMOD) for symmetric
# indefinite, then to LU (UMFPACK) as a last resort.
# Uses check=false to avoid expensive exception/stack-trace overhead.
#
# Note: diagonal perturbation (K + ε·diag) was tried but actually degrades
# solution quality on ill-conditioned systems.  LDLᵀ handles genuine
# near-singularity better than a perturbed Cholesky.
#
# After `update!(; values_only=true)` the old factor is kept but marked stale.
# If its sparsity-pattern key still matches, it is refactored in place with
# `cholesky!`, which skips CHOLMOD's symbolic analysis (fill-reducing ordering),
# the dominant cost for large shell models.
const _N_FACTORIZATIONS = Ref(0)

function _get_factorization(model::AbstractModel)
    old = model._factorization
    if old !== nothing && !model._factorization_stale
        return old
    end
    _N_FACTORIZATIONS[] += 1
    idx = model.freeDOFs
    K = Symmetric(model.S[idx, idx])
    key = _factor_pattern_key(model)
    if old isa SparseArrays.CHOLMOD.Factor && model._factor_key == key && key != 0
        cholesky!(old, K; check=false)
        if issuccess(old)
            model._factorization_stale = false
            # In-place numeric refactorization: reused symbolic analysis.
            log_analysis!(model, :refactor; n_free = length(idx))
            return old
        end
    end
    # Full factorization: new symbolic analysis (fill-reducing ordering).
    log_analysis!(model, :factor; n_free = length(idx))
    model._factorization = nothing
    model._factorization_stale = false
    model._factor_key = UInt(0)
    fact = cholesky(K; check=false)
    issuccess(fact) && (model._factor_key = key)
    if !issuccess(fact)
        @warn "Stiffness matrix not SPD — trying LDLᵀ"
        fact = ldlt(K; check=false)
        if !issuccess(fact)
            @warn "LDLᵀ failed — falling back to LU"
            K_unsym = model.S[idx, idx]
            fact = lu(K_unsym)
            # Spot-check: solve K*x = e₁ and verify residual isn't huge
            n = size(K_unsym, 1)
            if n > 0
                e1 = zeros(n); e1[1] = 1.0
                x1 = fact \ e1
                resid = norm(K_unsym * x1 - e1) / max(norm(e1), 1.0)
                if resid > 1e-6
                    @error "Stiffness matrix appears singular or severely ill-conditioned " *
                           "(residual = $(resid)). Check boundary conditions and supports."
                end
            end
        end
    end
    model._factorization = fact
    return fact
end

"""Key identifying the free-DOF sparsity pattern of `model.S` (same `S` object, nnz, and free DOFs)."""
_factor_pattern_key(model::AbstractModel) =
    hash((objectid(model.S), nnz(model.S), model.freeDOFs))

"""Mark the cached factorization stale after an in-place, same-pattern change to `S` values."""
function _invalidate_factor_values!(model::AbstractModel)
    model._factorization === nothing || (model._factorization_stale = true)
    return model
end

"""Key identifying the master/slave *reduced* sparsity pattern (`Krr`): the source
`S` object and nnz, the retained-free DOF set, and `Krr`'s own column/row structure."""
_reduced_factor_pattern_key(model::AbstractModel, Krr::Symmetric,
                            retained_free::AbstractVector{Int}) =
    hash((objectid(model.S), nnz(model.S), retained_free,
          Krr.data.colptr, Krr.data.rowval))

"""
    _get_reduced_factorization!(model, Krr::Symmetric, retained_free) -> factor

Factor the rigid-link master/slave *reduced* stiffness `Krr = Gᵀ S G` on the
retained free DOFs (Bathe §4.4), reusing the cached symbolic factorization when
the reduced sparsity pattern is unchanged — the same `cholesky!` trick
[`_get_factorization`](@ref) uses for the unreduced system, but applied to the
reduced matrix the rigid-link `solve!` builds.

The rigid-link solve path owns `model._factorization`/`_factorization_stale`/
`_factor_key` (it is only reached when links are present, where the unreduced
`_get_factorization` is never called), so the reduced factor is cached in those
slots and driven by the existing lifecycle: `update!(; values_only=true)` marks
it stale (→ in-place numeric refactor), a full `update!`/`process!` nulls it
(→ fresh factor). Every reuse is gated by [`_reduced_factor_pattern_key`](@ref);
any mismatch or non-SPD result falls back to a fresh `cholesky` → `ldlt` → `lu`,
identical to a cold factorization, so a reuse is only taken when it is
bit-identical to factoring from scratch.
"""
function _get_reduced_factorization!(model::AbstractModel, Krr::Symmetric,
                                     retained_free::AbstractVector{Int})
    key = _reduced_factor_pattern_key(model, Krr, retained_free)
    old = model._factorization
    # Same reduced pattern: in-place numeric refactorization from the CURRENT
    # `Krr` values, reusing the fill-reducing ordering from the symbolic
    # analysis. `cholesky!` always recomputes the numerics, so the returned
    # factor decomposes exactly this `Krr` no matter why it changed — the reuse
    # is of the (pattern-only) symbolic work, never of a stale numeric factor.
    if old isa SparseArrays.CHOLMOD.Factor && model._factor_key == key && key != 0
        cholesky!(old, Krr; check = false)
        if issuccess(old)
            model._factorization_stale = false
            log_analysis!(model, :refactor; n_free = length(retained_free))
            return old
        end
    end
    # Fresh factorization: new symbolic analysis, same SPD → LDLᵀ → LU fallback
    # ladder as the unreduced path and the previous rigid-link code.
    _N_FACTORIZATIONS[] += 1
    log_analysis!(model, :factor; n_free = length(retained_free))
    model._factorization = nothing
    model._factorization_stale = false
    model._factor_key = UInt(0)
    fact = cholesky(Krr; check = false)
    if issuccess(fact)
        model._factor_key = key
    else
        @warn "RigidLink reduced stiffness not SPD — trying LDLᵀ"
        fact = ldlt(Krr; check = false)
        if !issuccess(fact)
            @warn "LDLᵀ failed — falling back to LU"
            fact = lu(Krr.data)
        end
    end
    model._factorization = fact
    return fact
end

# =============================================================================
# #7 phase 2b — Caller-owned reduced-factorization cache (lateral solves)
# =============================================================================

"""
    RigidLinkFactorCache()

Caller-owned cache for the reduced rigid-link factorization used by the
non-mutating [`solve`](@ref)`(model, cases, links; factor_cache=…)`.

The mutating `solve!`/plain-`solve` reuse path caches a *single* factor in
`model._factorization`, which cannot serve the `geometric_stiffness` (P-Δ) or
`exclude_shells` (frame-only drift) variants: each builds a *different* reduced
matrix `Krr`, and the lateral stage interleaves them across many load cases and
gravity combinations. This cache holds one CHOLMOD symbolic factorization per
distinct reduced pattern (keyed by [`_reduced_factor_pattern_key`](@ref) plus
the solve variant), so a sweep that changes only *values* — the story geometric
stiffness `K_G` across gravity keys, or a section update between sizing
iterations — reuses the fill-reducing ordering via `cholesky!` instead of a
fresh symbolic analysis every call.

Pass it as the `factor_cache` keyword; omit it (or pass `nothing`) for the
previous fresh-every-call behaviour. The key folds in `objectid(model.S)` and
`nnz(model.S)`, so a remesh (new `S`) never matches a stale entry; still, the
owner should discard or clear the cache when the mesh/topology changes to avoid
holding dead factors. Only SPD Cholesky factors are cached — an indefinite
`Krr` (e.g. a buckling P-Δ system) falls back to LDLᵀ/LU every call, exactly as
the uncached path does.
"""
mutable struct RigidLinkFactorCache
    factors::Dict{UInt, SparseArrays.CHOLMOD.Factor}
end
RigidLinkFactorCache() = RigidLinkFactorCache(Dict{UInt, SparseArrays.CHOLMOD.Factor}())

"""Drop every cached factor (call when the mesh/topology changes)."""
clear!(cache::RigidLinkFactorCache) = (empty!(cache.factors); cache)

"""
    _cached_reduced_factor!(cache, model, Krr::Symmetric, retained_free, key) -> factor

Factor the reduced stiffness `Krr` using `cache`, reusing the cached CHOLMOD
symbolic analysis via `cholesky!` when an entry for `key` already exists. Like
[`_get_reduced_factorization!`](@ref) this always recomputes the numerics from
the *current* `Krr`, so a reuse is bit-identical to factoring from scratch; only
the (pattern-only) fill-reducing ordering is reused. `key` must encode both the
reduced sparsity pattern and the solve variant so two matrices with coincident
patterns but different semantics (e.g. `S` vs `S + K_G`) never share a factor.
Only SPD factors are cached; an indefinite `Krr` falls back to LDLᵀ → LU and is
not stored, so the cache never holds a factor it cannot refresh in place.
"""
function _cached_reduced_factor!(cache::RigidLinkFactorCache, model::AbstractModel,
                                 Krr::Symmetric, retained_free::AbstractVector{Int},
                                 key::UInt)
    old = get(cache.factors, key, nothing)
    # Same reduced pattern + variant: in-place numeric refactorization from the
    # CURRENT `Krr`, reusing the cached fill-reducing ordering.
    if old !== nothing && key != 0
        cholesky!(old, Krr; check = false)
        if issuccess(old)
            log_analysis!(model, :refactor; n_free = length(retained_free))
            return old
        end
    end
    # Fresh symbolic analysis. Same SPD → LDLᵀ → LU ladder as the uncached path.
    # Like the rest of the non-mutating `solve`, this is invisible to
    # `_N_FACTORIZATIONS`; the `:factor`/`:refactor` log events are the reuse
    # signal.
    log_analysis!(model, :factor; n_free = length(retained_free))
    fact = cholesky(Krr; check = false)
    if issuccess(fact)
        key != 0 && (cache.factors[key] = fact)
        return fact
    end
    # Indefinite reduced stiffness (e.g. a buckling P-Δ system): do not cache a
    # factor we cannot `cholesky!`-refresh. Drop any stale SPD entry for `key`.
    delete!(cache.factors, key)
    @warn "RigidLink reduced stiffness not SPD — trying LDLᵀ"
    fact = ldlt(Krr; check = false)
    if !issuccess(fact)
        @warn "LDLᵀ failed — falling back to LU"
        fact = lu(Krr.data)
    end
    return fact
end

