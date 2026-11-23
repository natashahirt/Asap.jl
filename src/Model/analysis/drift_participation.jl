# =============================================================================
# Drift participation — virtual-work (unit-load) sensitivity
# =============================================================================
#
# For a solved model, the lateral displacement δ at a chosen free DOF
# decomposes into per-element contributions by the unit-load theorem. Apply a
# unit virtual load at that DOF; reusing the cached stiffness factorization
# (the virtual load is just one more right-hand side — no refactorization), the
# virtual displacement field u_v gives each frame element's participation
#
#     c_e = u_{v,e}ᵀ K_e u_{r,e}           (K_e the element's GCS stiffness)
#
# and, by reciprocity (Maxwell–Betti),
#
#     Σ_e c_e = u_vᵀ S u_r = u_vᵀ F_r = u_r[dof] = δ.
#
# This is exact and free of any moment sign convention: it reuses only the
# assembled element stiffness `element.K`. `c_e` is the element's share of the
# drift and is the drift-participation coefficient the optimality-criteria
# resizer consumes (δ = Σ c_e, and ∂δ/∂(stiffness)_e = −c_e / stiffness_e, so
# material goes where participation-per-unit-weight is highest).
#
# For a frame-only model `Σ c_e == δ` to rounding (`nonframe ≈ 0`), which is a
# built-in correctness check. For a model whose core/walls are shells (a dual
# system), `nonframe = δ − Σ c_e` is the non-frame (core) share of the drift —
# a meaningful quantity, not an error.
#
# References: Maxwell–Mohr unit-load theorem; Baker (1990), energy-based design
# of lateral systems; Chan, Grierson & Sherbourne (1995), optimality criteria
# for drift-controlled tall steel buildings.

"""
    DriftParticipation

Result of [`drift_participation`](@ref):

- `dof`      — the queried free global DOF.
- `δ`        — the drift (real displacement) at `dof`.
- `element_ids` — frame-element IDs, aligned with `c`.
- `c`        — per-element drift-participation coefficients (`Σ c ≈ δ` for a
               frame-only model).
- `nonframe` — `δ − Σ c`: ~0 for a frame-only model (a correctness check), or
               the non-frame (shell/core, spring, rigid-link) share of the
               drift for a mixed model.
"""
struct DriftParticipation
    dof::Int
    δ::Float64
    element_ids::Vector{Int64}
    c::Vector{Float64}
    nonframe::Float64
end

# Frame elements that carry a 12-DOF `globalID` and GCS `K`.
_participation_elements(model::FrameModel) = model.elements
_participation_elements(model::TrussModel) = model.elements
_participation_elements(model::Model)      = model.frame_elements

"""
    drift_participation(model, dof::Integer) -> DriftParticipation

Decompose the lateral displacement at free global DOF `dof` into per-element
drift-participation coefficients `c_e = u_vᵀ K_e u_r` by the unit-load theorem,
reusing the model's cached stiffness factorization (the virtual unit load is one
extra right-hand side — no refactorization).

`model` must already be solved (`solve!`, so `model.u` holds the real
displacement field) and `dof` must be a free DOF. `Σ c_e == δ == model.u[dof]`
by reciprocity; the shortfall is reported as `nonframe` (≈0 for a frame-only
model, else the core/shell share).

    drift_participation(model, node::Node, direction::Integer)

Convenience overload: `direction` is the node's local DOF (1–3 translations,
4–6 rotations), e.g. `1` for global-X story drift.
"""
function drift_participation(model::AbstractModel, dof::Integer)
    d = Int(dof)
    idx = model.freeDOFs
    pos = findfirst(==(d), idx)
    pos === nothing && error(
        "drift_participation: DOF $d is not a free DOF (it may be restrained " *
        "or out of range). Query a free lateral DOF.")
    length(model.u) == model.nDOFs || error(
        "drift_participation: model is not solved; call solve!(model) first.")

    # Virtual unit load at `dof`, solved on the cached factor (extra RHS only).
    fact = _get_factorization(model)
    ev = zeros(length(idx))
    ev[pos] = 1.0
    Uv = fact \ ev
    uv = zeros(model.nDOFs)
    uv[idx] = Uv
    ur = model.u

    els = _participation_elements(model)
    ids = Vector{Int64}(undef, length(els))
    c = Vector{Float64}(undef, length(els))
    @inbounds for (i, el) in enumerate(els)
        g = el.globalID
        c[i] = dot(view(uv, g), el.K * view(ur, g))
        ids[i] = el.elementID
    end

    δ = ur[d]
    nonframe = δ - sum(c)
    return DriftParticipation(d, δ, ids, c, nonframe)
end

function drift_participation(model::AbstractModel, node::Node, direction::Integer)
    1 <= direction <= 6 || error(
        "drift_participation: direction must be 1–6 (1–3 translation, 4–6 " *
        "rotation); got $direction.")
    return drift_participation(model, node.globalID[direction])
end
