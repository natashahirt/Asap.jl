# =============================================================================
# Drift participation — virtual-work (unit-load) sensitivity
# =============================================================================
#
# For a solved model, a lateral displacement functional δ = vᵀu_r (the drift at
# a DOF, an inter-story drift, or any linear combination of DOFs selected by the
# virtual-load vector `v`) decomposes into per-element contributions by the
# unit-load theorem. Apply the virtual load `v`; reusing the cached stiffness
# factorization (`v` is just one more right-hand side — no refactorization), the
# virtual displacement field u_v gives each frame element's participation
#
#     c_e = u_{v,e}ᵀ K_e u_{r,e}           (K_e the element's GCS stiffness)
#
# and, by reciprocity (Maxwell–Betti),
#
#     Σ_e c_e = u_vᵀ S u_r = vᵀ u_r = δ.
#
# This is exact and free of any moment sign convention: it reuses only the
# assembled element stiffness `element.K`. `c_e` is the element's share of the
# drift functional and is the drift-participation coefficient the
# optimality-criteria resizer consumes (δ = Σ c_e, and ∂δ/∂(stiffness)_e =
# −c_e / stiffness_e, so material goes where participation-per-unit-weight is
# highest).
#
# For a frame-only model `Σ c_e == δ` to rounding (`nonframe ≈ 0`), a built-in
# correctness check. For a model whose core/walls are shells (a dual system),
# `nonframe = δ − Σ c_e` is the non-frame (core) share of the drift — a
# meaningful quantity, not an error, and exactly the "core background" the dual
# drift governor sizes the frame against (δ_eff = δ_allow − nonframe).
#
# Virtual loads supported:
#   * a single free DOF            → absolute displacement at that DOF;
#   * a differential DOF pair      → inter-story drift u[i] − u[i−1] (the
#                                    governing ASCE 7-22 §12.8.6 limit);
#   * an arbitrary vector `v`      → any linear displacement functional;
#   * a matrix `V` (columns = v_k) → several functionals from ONE multi-RHS
#                                    back-substitution (per-story multi-λ sizing).
#
# Rigid-link (rigid-diaphragm) models are solved in the reduced master/slave
# space: the virtual load is mapped through Gᵀ, solved on the cached *reduced*
# factor, and the virtual field recovered through G — the same condensation
# `solve!(model, links)` uses — so the diaphragm-master story DOFs participate
# correctly. Pass the model's `links` for this path.
#
# References: Maxwell–Mohr unit-load theorem; Baker (1990), energy-based design
# of lateral systems; Chan, Grierson & Sherbourne (1995), optimality criteria
# for drift-controlled tall steel buildings.

"""
    DriftParticipation

Result of [`drift_participation`](@ref) for one virtual load:

- `δ`        — the drift functional value `vᵀ u_r` (real displacement).
- `element_ids` — frame-element IDs, aligned with `c`.
- `c`        — per-element drift-participation coefficients (`Σ c ≈ δ` for a
               frame-only model).
- `nonframe` — `δ − Σ c`: ~0 for a frame-only model (a correctness check), or
               the non-frame (shell/core, spring, rigid-link) share of the
               drift for a mixed model.
"""
struct DriftParticipation
    δ::Float64
    element_ids::Vector{Int64}
    c::Vector{Float64}
    nonframe::Float64
end

"""
    BatchedDriftParticipation

Result of [`drift_participation`](@ref) for several virtual loads (a matrix `V`
whose columns are the virtual-load vectors `v_k`), from one multi-RHS solve:

- `δ`        — drift functional value for each column, `δ[k] = v_kᵀ u_r`.
- `element_ids` — frame-element IDs, aligned with the rows of `c`.
- `c`        — `c[e, k]`, element `e`'s participation in functional `k`.
- `nonframe` — `nonframe[k] = δ[k] − Σ_e c[e, k]`, the non-frame share of each.
"""
struct BatchedDriftParticipation
    δ::Vector{Float64}
    element_ids::Vector{Int64}
    c::Matrix{Float64}
    nonframe::Vector{Float64}
end

# Frame elements that carry a 12-DOF `globalID` and GCS `K`.
_participation_elements(model::FrameModel) = model.elements
_participation_elements(model::TrussModel) = model.elements
_participation_elements(model::Model)      = model.frame_elements

"""Rigid-link reduction `(G, retained_free)` for `links`, or `nothing` when the
links tie no DOFs (the unreduced factor is then the right one)."""
function _participation_reduction(model::AbstractModel,
                                  links::Union{Nothing, AbstractVector})
    (links === nothing || isempty(links)) && return nothing
    G, dependent = _build_rigidlink_reduction(model, links)
    isempty(dependent) && return nothing
    retained_free = setdiff(model.freeDOFs, dependent)
    return (G = G, retained_free = retained_free)
end

"""
    _virtual_fields(model, V, reduction) -> Matrix

Virtual displacement fields (full `nDOFs × k`) for the virtual-load columns of
`V` (full `nDOFs × k`), solved on the cached factorization. `reduction` is the
rigid-link `(G, retained_free)` or `nothing` for the unreduced system. One
multi-RHS back-substitution covers every column.
"""
function _virtual_fields(model::AbstractModel, V::AbstractMatrix,
                         reduction::Union{Nothing, NamedTuple})
    uv = zeros(model.nDOFs, size(V, 2))
    if reduction === nothing
        idx = model.freeDOFs
        fact = _get_factorization(model)
        Uv = fact \ Matrix(V[idx, :])
        uv[idx, :] = Uv
    else
        G = reduction.G
        rf = reduction.retained_free
        # Reduce the virtual load (Gᵀ v), solve on the retained free DOFs, and
        # recover the full virtual field u_v = G u_retained — the condensation
        # `solve!(model, links)` uses (Bathe §4.4).
        S = model.S
        Krr = Symmetric((transpose(G) * S * G)[rf, rf])
        fact = _get_reduced_factorization!(model, Krr, rf)
        Fc = transpose(G) * V
        Ur = fact \ Matrix(Fc[rf, :])
        full = zeros(model.nDOFs, size(V, 2))
        full[rf, :] = Ur
        uv = G * full
    end
    return uv
end

"""
    _participation(model, V, δ, reduction) -> (ids, c, nonframe)

Per-element participation `c[e, k] = u_{v,e}(k)ᵀ K_e u_{r,e}` for every virtual
field (column `k`) and the non-frame remainder `δ[k] − Σ_e c[e, k]`.
"""
function _participation(model::AbstractModel, V::AbstractMatrix,
                        δ::AbstractVector, reduction::Union{Nothing, NamedTuple})
    uv = _virtual_fields(model, V, reduction)
    ur = model.u
    els = _participation_elements(model)
    nk = size(V, 2)
    ids = Vector{Int64}(undef, length(els))
    c = Matrix{Float64}(undef, length(els), nk)
    @inbounds for (i, el) in enumerate(els)
        g = el.globalID
        Kur = el.K * view(ur, g)              # real element end forces, once
        for k in 1:nk
            c[i, k] = dot(view(uv, g, k), Kur)
        end
        ids[i] = el.elementID
    end
    nonframe = [δ[k] - sum(view(c, :, k)) for k in 1:nk]
    return ids, c, nonframe
end

"""Column `e_dof` (unit load at one DOF) as a sparse `nDOFs`-vector."""
function _unit_load(model::AbstractModel, dof::Integer)
    v = zeros(model.nDOFs)
    v[dof] = 1.0
    return v
end

"""
    _require_free(model, dof)

Throw when `dof` is not a free DOF (restrained or out of range), so a virtual
load is applied only where the structure can displace.
"""
function _require_free(model::AbstractModel, dof::Integer)
    (1 <= dof <= model.nDOFs && dof in model.freeDOFs) || error(
        "drift_participation: DOF $dof is not a free DOF (it may be restrained " *
        "or out of range). Query a free lateral DOF.")
    return nothing
end

"""
    drift_participation(model, v::AbstractVector; links=nothing) -> DriftParticipation

Decompose the displacement functional `δ = vᵀ u_r` into per-element
participation `c_e = u_vᵀ K_e u_r` by the unit-load theorem, reusing the model's
cached factorization (the virtual load `v`, an `nDOFs`-vector, is one extra
right-hand side — no refactorization). `model` must already be solved.

Pass `links` (the model's rigid links) for a rigid-diaphragm model: the virtual
load is then solved on the reduced master/slave system, matching how the model
was solved. `Σ c_e == δ` for a frame-only model; the shortfall is `nonframe`
(the core/shell share for a dual system).
"""
function drift_participation(model::AbstractModel, v::AbstractVector{<:Real};
                             links::Union{Nothing, AbstractVector} = nothing)
    length(v) == model.nDOFs || error(
        "drift_participation: virtual load has length $(length(v)); expected " *
        "nDOFs = $(model.nDOFs).")
    length(model.u) == model.nDOFs || error(
        "drift_participation: model is not solved; call solve!(model) first.")
    reduction = _participation_reduction(model, links)
    V = reshape(collect(Float64, v), model.nDOFs, 1)
    δ = [dot(v, model.u)]
    ids, c, nonframe = _participation(model, V, δ, reduction)
    return DriftParticipation(δ[1], ids, vec(c), nonframe[1])
end

"""
    drift_participation(model, dof::Integer; links=nothing) -> DriftParticipation

Participation of the absolute displacement at free DOF `dof` (virtual load = a
unit load at `dof`).
"""
function drift_participation(model::AbstractModel, dof::Integer;
                             links::Union{Nothing, AbstractVector} = nothing)
    _require_free(model, Int(dof))
    return drift_participation(model, _unit_load(model, Int(dof)); links = links)
end

"""
    drift_participation(model, dof_top, dof_bottom; links=nothing) -> DriftParticipation

Participation of the **inter-story drift** `u[dof_top] − u[dof_bottom]` (virtual
load = +1 at `dof_top`, −1 at `dof_bottom`), the governing ASCE 7-22 §12.8.6
quantity. `dof_bottom` may be a restrained base DOF (its virtual term then does
no work); `dof_top` must be free.
"""
function drift_participation(model::AbstractModel, dof_top::Integer, dof_bottom::Integer;
                             links::Union{Nothing, AbstractVector} = nothing)
    _require_free(model, Int(dof_top))
    v = _unit_load(model, Int(dof_top))
    1 <= Int(dof_bottom) <= model.nDOFs || error(
        "drift_participation: dof_bottom $(dof_bottom) is out of range.")
    v[Int(dof_bottom)] -= 1.0
    return drift_participation(model, v; links = links)
end

"""
    drift_participation(model, node::Node, direction::Integer; links=nothing)

Convenience overload: `direction` is the node's local DOF (1–3 translations,
4–6 rotations), e.g. `1` for global-X story drift at a diaphragm master.
"""
function drift_participation(model::AbstractModel, node::Node, direction::Integer;
                             links::Union{Nothing, AbstractVector} = nothing)
    1 <= direction <= 6 || error(
        "drift_participation: direction must be 1–6 (1–3 translation, 4–6 " *
        "rotation); got $direction.")
    return drift_participation(model, node.globalID[direction]; links = links)
end

"""
    drift_participation(model, V::AbstractMatrix; links=nothing) -> BatchedDriftParticipation

Participation for several virtual loads at once: the columns of `V` (each an
`nDOFs`-vector) are solved in ONE multi-RHS back-substitution, returning
`c[e, k]` for element `e` and functional `k`. This is the per-story multi-λ
sizing path — 40–80 story/direction drift constraints at the cost of one solve
rather than one per constraint.
"""
function drift_participation(model::AbstractModel, V::AbstractMatrix{<:Real};
                             links::Union{Nothing, AbstractVector} = nothing)
    size(V, 1) == model.nDOFs || error(
        "drift_participation: virtual-load matrix has $(size(V, 1)) rows; " *
        "expected nDOFs = $(model.nDOFs).")
    length(model.u) == model.nDOFs || error(
        "drift_participation: model is not solved; call solve!(model) first.")
    reduction = _participation_reduction(model, links)
    Vf = Matrix{Float64}(V)
    δ = [dot(view(Vf, :, k), model.u) for k in 1:size(Vf, 2)]
    ids, c, nonframe = _participation(model, Vf, δ, reduction)
    return BatchedDriftParticipation(δ, ids, c, nonframe)
end

# =============================================================================
# Linear response — reuse the cached factor for an arbitrary load (first-order
# sizing sensitivity)
# =============================================================================

"""
    linear_response(model, f::AbstractVector; links=nothing) -> Vector

Full displacement response `u = K⁻¹ f` to a nodal load vector `f` (length
`nDOFs`), reusing the model's cached stiffness factorization — one
back-substitution, no refactorization. For a rigid-diaphragm model pass `links`;
`f` is then reduced through `Gᵀ`, solved on the cached reduced factor, and the
response recovered through `G` (`u = G u_retained`). `model` must already be
solved/factored.

This is the primitive behind first-order drift-sizing sensitivity: under a fixed
external load, perturbing the element stiffnesses by `ΔK` shifts the
displacement by

    Δu ≈ −K⁻¹ (ΔK · u) = linear_response(model, −ΔK*u; links),

so a resize step's redistributed displacements (and hence member demands) can be
predicted from the *current* factor without a fresh solve — one extra RHS, like
the virtual-load solve in [`drift_participation`](@ref). Sanity check: feeding
the real applied load reproduces the real displacement,
`linear_response(model, model.P .- model.Pf) ≈ model.u`.
"""
function linear_response(model::AbstractModel, f::AbstractVector{<:Real};
                         links::Union{Nothing, AbstractVector} = nothing)
    length(f) == model.nDOFs || error(
        "linear_response: load vector has length $(length(f)); expected " *
        "nDOFs = $(model.nDOFs).")
    length(model.u) == model.nDOFs || error(
        "linear_response: model is not solved; call solve!(model) first.")
    reduction = _participation_reduction(model, links)
    F = reshape(collect(Float64, f), model.nDOFs, 1)
    return vec(_virtual_fields(model, F, reduction))
end
