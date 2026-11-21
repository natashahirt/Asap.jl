#=
Rigid Link Element for Asap.jl
==============================

`RigidLink` is a two-node kinematic-offset primitive: it ties a *slave* node
(`nodeB`) to a *master* node (`nodeA`) through a rigid body so that the slave
follows the master's motion exactly (up to selectively released DOFs). It is
NOT a stiffness-bearing element — it contributes no `K`. Instead it is resolved
during the solve by *master/slave DOF reduction*: the constrained slave DOFs are
algebraically substituted out of the global system before the solve, and their
displacements are recovered by back-substitution afterwards.

Kinematics (small-rotation rigid-body transformation):
    u_B = u_A + θ_A × r_AB
    θ_B = θ_A
with r_AB = position_B − position_A (the offset from master to slave).

Why reduction (not penalty springs or Lagrange multipliers):
- A stiff penalty spring introduces a large stiffness ratio that ill-conditions
  the global matrix; the answer depends on an arbitrary penalty value.
- Lagrange multipliers enlarge the system and make it indefinite (no longer SPD),
  losing the fast Cholesky path.
- Master/slave reduction is exact, keeps the reduced system SPD, and shrinks it.
  See Bathe, "Finite Element Procedures" §4.4 (constraint / rigid-link handling)
  and Cook et al., "Concepts and Applications of Finite Element Analysis" §6.3.

Applications: rigid diaphragm master nodes, eccentric brace/gusset working
points, composite-slab kinematic coupling, coupled-wall link beams, and hinged
connections at offsets (via the per-DOF release vector).
=#

using LinearAlgebra: I, dot, Symmetric, cholesky, ldlt, lu, issuccess
using SparseArrays: sparse, SparseMatrixCSC, rowvals, nonzeros

"""
    RigidLink(nodeA, nodeB; release=(false,false,false,false,false,false), id=:rigidlink)

Create a rigid kinematic link tying slave `nodeB` to master `nodeA`.

# Arguments
- `nodeA::Node`: master node — the retained node that drives the motion.
- `nodeB::Node`: slave node — its DOFs are constrained to follow `nodeA`.

# Keyword Arguments
- `release::NTuple{6, Bool}`: per-DOF relaxation of the constraint, ordered
  `(u_x, u_y, u_z, θ_x, θ_y, θ_z)` in the global coordinate system.
  `false` = rigid (tied to master), `true` = released (that slave DOF stays a
  free, independent DOF).
- `id::Symbol`: optional identifier (default `:rigidlink`).

# Release examples
- `(false, false, false, false, false, false)` — fully rigid link.
- `(false, false, false, false, true, true)` — pin-and-collar: rigid in
  translation, free to rotate about the global y and z axes.
- `(false, false, false, true, false, false)` — torsion release: rigid
  translation and y/z rotation, free about the global x axis.

# Notes
- The slave node must be a regular `Node` present in the model's node list; its
  fixities are ignored on any tied DOF (the tie governs), and a warning is
  emitted if a tied slave DOF is also a support.
- Chained links (a master node that is itself another link's slave) are not
  supported and raise an error.
"""
mutable struct RigidLink <: AbstractElement
    nodeA::Node                     # master
    nodeB::Node                     # slave
    release::NTuple{6, Bool}        # per-DOF: true = released (free), false = rigid
    elementID::Int64
    id::Symbol
end

function RigidLink(nodeA::Node, nodeB::Node;
                   release::NTuple{6, Bool} = (false, false, false, false, false, false),
                   id::Symbol = :rigidlink)
    return RigidLink(nodeA, nodeB, release, 0, id)
end

"""
    rigid_link_transformation(nodeA, nodeB) -> Matrix{Float64}

Return the 6×6 rigid-body transformation `T` mapping the master node's DOFs to
the slave node's DOFs in the global coordinate system: `[u_B; θ_B] = T [u_A; θ_A]`.

Derived from the standard 6×6 rigid-offset (Te) matrix used in FE textbooks
(Bathe §4.4; Cook §6.3) and mirrors `eccentricity_transformation!` in `K.jl`.

Small-rotation rigid-body transformation:  u_B = u_A + θ_A × r_AB,  θ_B = θ_A,
with r_AB = (rx, ry, rz) = position_B − position_A. Expanding θ_A × r_AB:
    u_Bx = u_Ax + θ_Ay·rz − θ_Az·ry
    u_By = u_Ay + θ_Az·rx − θ_Ax·rz
    u_Bz = u_Az + θ_Ax·ry − θ_Ay·rx
"""
function rigid_link_transformation(nodeA::Node, nodeB::Node)
    rx = to_meters(nodeB.position[1]) - to_meters(nodeA.position[1])
    ry = to_meters(nodeB.position[2]) - to_meters(nodeA.position[2])
    rz = to_meters(nodeB.position[3]) - to_meters(nodeA.position[3])

    T = Matrix{Float64}(I, 6, 6)
    # Translation ← master rotation coupling (θ_A × r_AB); rotation rows are identity.
    T[1, 5] = rz;   T[1, 6] = -ry
    T[2, 4] = -rz;  T[2, 6] = rx
    T[3, 4] = ry;   T[3, 5] = -rx
    return T
end

"""
    _build_rigidlink_reduction(model, links) -> (G, dependent)

Build the master/slave reduction transformation `G` (nDOFs × nDOFs) and the list
of `dependent` (slave-constrained) global DOF indices.

`G` maps the reduced (retained-DOF) displacement vector into the full DOF space:
`u_full = G · u_retained`. It is the identity on retained DOFs; each dependent
(tied slave) DOF row holds the master-DOF combination from
`rigid_link_transformation`, so `u_slave = T · u_master` is enforced exactly.
"""
function _build_rigidlink_reduction(model::AbstractModel, links::AbstractVector{RigidLink})
    n = model.nDOFs

    # combos[d] = list of (master_dof, coefficient) for each dependent (tied) slave DOF d
    combos = Dict{Int, Vector{Tuple{Int, Float64}}}()
    dependent = Int[]

    for link in links
        T = rigid_link_transformation(link.nodeA, link.nodeB)
        gA = link.nodeA.globalID
        gB = link.nodeB.globalID
        for i in 1:6
            link.release[i] && continue          # released: slave DOF stays independent
            d = gB[i]
            if !model.DOFs[d]
                @warn "RigidLink: slave DOF $d is a support (fixed); leaving it as a support and skipping the rigid tie for that DOF."
                continue
            end
            haskey(combos, d) && error("RigidLink: slave DOF $d is constrained by more than one link (overconstrained).")
            entries = Tuple{Int, Float64}[]
            for j in 1:6
                c = T[i, j]
                c == 0.0 && continue
                push!(entries, (gA[j], c))
            end
            combos[d] = entries
            push!(dependent, d)
        end
    end

    depset = Set(dependent)

    # Assemble G in COO form: identity on retained rows, master-combination on dependent rows.
    Is = Int[]; Js = Int[]; Vs = Float64[]
    for k in 1:n
        if haskey(combos, k)
            for (m, c) in combos[k]
                # Masters must themselves be retained; chained links are not supported.
                m in depset && error("RigidLink: chained links are not supported (master DOF $m is itself a tied slave).")
                push!(Is, k); push!(Js, m); push!(Vs, c)
            end
        else
            push!(Is, k); push!(Js, k); push!(Vs, 1.0)
        end
    end

    G = sparse(Is, Js, Vs, n, n)
    return G, dependent
end

"""
    solve!(model, links::AbstractVector{RigidLink}; reprocess=true, postprocess=:all)

Solve a frame/unified model subject to rigid kinematic links via master/slave
DOF reduction (Bathe §4.4; Cook §6.3).

Steps:
1. `process!` the model to assemble the full stiffness `S` and load vectors.
2. Build the reduction `G` from the links (see `_build_rigidlink_reduction`).
3. Condense: `S_c = Gᵀ S G`, `F_c = Gᵀ (P − Pf)`, then solve the reduced system
   on the retained free DOFs (all free DOFs minus the tied slave DOFs).
4. Back-substitute the slave displacements via `u = G · u_retained`.
5. Post-process using the *full* `S` and the recovered *full* `u`, so element
   forces and reactions reflect the true physical (constrained) state.

Only free slave DOFs are tied; released DOFs (`release[i] == true`) remain as
ordinary free DOFs in the reduced system.
"""
function solve!(model::AbstractModel, links::AbstractVector{RigidLink};
                reprocess::Bool = true, postprocess::Symbol = :all)
    if !model.processed || reprocess
        process!(model)
    end

    if isempty(links)
        solve!(model; postprocess = postprocess)
        return model
    end

    G, dependent = _build_rigidlink_reduction(model, links)

    if isempty(dependent)
        solve!(model; postprocess = postprocess)
        return model
    end

    # Full physical system (kept intact for post-processing).
    S = model.S
    F = model.P .- model.Pf

    # Master/slave reduction via kinematic transformation (Bathe §4.4, Cook §6.3):
    # substitute slave DOFs out of the global system before the solve.
    Gt = transpose(G)
    Sc = Gt * S * G
    Fc = Gt * F

    retained_free = setdiff(model.freeDOFs, dependent)

    Krr = Symmetric(Sc[retained_free, retained_free])
    # #7: factor the reduced system through the model's factor cache. On a
    # value-only re-solve (same mesh/links/geometry; only section stiffness
    # changed) `Krr`'s sparsity pattern is unchanged, so `cholesky!` reuses the
    # symbolic analysis instead of a fresh factorization. The cache is exact
    # (pattern-keyed, SPD-guarded) and falls back to a fresh factorization
    # otherwise — see `_get_reduced_factorization!`.
    fact = _get_reduced_factorization!(model, Krr, retained_free)

    Ur = fact \ Fc[retained_free]

    u = zeros(model.nDOFs)
    u[retained_free] = Ur

    # Back-substitute slave displacements: u_slave = T · u_master (i.e. u = G·u).
    model.u = G * u

    model.compliance = dot(Ur, Fc[retained_free])
    # Keep the reduced factor cached in `model._factorization` for reuse by the
    # next value-only re-solve; the lifecycle (`update!`/`process!`) invalidates
    # it when `S` changes.

    post_process!(model; targets = postprocess)
    return model
end

"""Degrees of freedom that belong to a shell node and to no frame element."""
function _shell_only_dofs(model::Model)
    frame_ids = Set{UInt}()
    for el in model.frame_elements
        push!(frame_ids, objectid(el.nodeStart), objectid(el.nodeEnd))
    end
    dofs = Int[]
    seen = Set{UInt}()
    for el in model.shell_elements
        for nd in el.nodes
            id = objectid(nd)
            (id in frame_ids || id in seen) && continue
            push!(seen, id)
            append!(dofs, nd.globalID)
        end
    end
    return dofs
end

"""
    solve(model, cases, links; geometric_stiffness=nothing, exclude_shells=false) -> Vector{Vector{Float64}}

Non-mutating multi-case solve on the rigid-link reduced system.

Same condensation as [`solve!(model, links)`](@ref): `G` from
`_build_rigidlink_reduction`, then `K_c = Gᵀ K G` on the retained free DOFs.
One factorization, one column of `F` per case. Each result is a full
displacement vector, slaves included (`u = G u_retained`), in metres.

Does not write `model.u`, `model.P`, or `model._factorization`, and does not
post-process. An empty link set, or a link set that ties nothing, uses
[`solve(model, cases)`](@ref) on the unreduced stiffness.

`geometric_stiffness`, when given, is an `nDOFs × nDOFs` symmetric matrix
added to the elastic stiffness before the reduction: the solve is on
`Gᵀ (S + K_G) G`. It is signed, so a story in compression passes a negative
`K_G` (Wilson & Habibullah 1987, Eq. 2: `K* = K − L`). It must not touch a
tied slave DOF; the caller puts story terms on the masters. The returned
displacements then satisfy equilibrium in the deformed position, and element
forces recovered as `K_e u_e` carry the second-order actions.

`exclude_shells` subtracts the shell stiffness and holds every degree of
freedom that only a shell touches. Diaphragm slaves stay tied. The numbering
does not change.
"""
function solve(model::Model,
               cases::Vector{<:Vector{<:AbstractLoad}},
               links::AbstractVector{RigidLink};
               geometric_stiffness::Union{Nothing, AbstractMatrix} = nothing,
               exclude_shells::Bool = false)
    model.processed || process!(model)
    ncases = length(cases)
    ncases == 0 && return Vector{Float64}[]

    normalized = Vector{Vector{AbstractLoad}}(undef, ncases)
    for j in 1:ncases
        normalized[j] = Vector{AbstractLoad}(cases[j])
    end

    if !exclude_shells && geometric_stiffness === nothing && isempty(links)
        return solve(model, normalized)
    end

    G, dependent = _build_rigidlink_reduction(model, links)
    if !exclude_shells && geometric_stiffness === nothing && isempty(dependent)
        return solve(model, normalized)
    end

    S = model.S
    # Drop the shell stiffness and hold shell-only nodes. The node numbering
    # stays the one process! already wrote, so a frame-only displacement lines
    # up with a displacement of the full model. A shell node tied to a rigid
    # diaphragm stays a slave: it adds no stiffness once the shells are gone.
    if exclude_shells && !isempty(model.shell_elements)
        S = S - assemble_stiffness(model.shell_elements, model.nDOFs)
    end
    held = exclude_shells ? setdiff(_shell_only_dofs(model), dependent) : Int[]
    if geometric_stiffness !== nothing
        size(geometric_stiffness) == size(S) || throw(DimensionMismatch(
            "geometric_stiffness is $(size(geometric_stiffness)); the model has $(size(S))."))
        Kg = sparse(geometric_stiffness)
        depset = Set(dependent)
        rows = rowvals(Kg)
        vals = nonzeros(Kg)
        for k in eachindex(rows)
            (vals[k] != 0 && rows[k] in depset) && throw(ArgumentError(
                "geometric_stiffness has terms on tied slave DOF $(rows[k]). " *
                "Put story terms on the master."))
        end
        S = S + Kg
    end
    Gt = transpose(G)
    Sc = Gt * S * G
    retained_free = setdiff(model.freeDOFs, dependent, held)
    Krr = Symmetric(Sc[retained_free, retained_free])
    # Rigid-diaphragm multi-case solve: one fresh reduced factorization reused
    # across every load case. Invisible to _N_FACTORIZATIONS.
    log_analysis!(model, :factor; n_free = length(retained_free), n_rhs = ncases)
    fact = cholesky(Krr; check = false)
    if !issuccess(fact)
        @warn "RigidLink reduced stiffness not SPD — trying LDLᵀ"
        fact = ldlt(Krr; check = false)
        if !issuccess(fact)
            @warn "LDLᵀ failed — falling back to LU"
            fact = lu(Sc[retained_free, retained_free])
        end
    end

    nret = length(retained_free)
    Fr = Matrix{Float64}(undef, nret, ncases)
    for j in 1:ncases
        Fj = create_F(model, normalized[j])
        Fr[:, j] = (Gt * Fj)[retained_free]
    end
    Ur = fact \ Fr

    n = model.nDOFs
    results = Vector{Vector{Float64}}(undef, ncases)
    for j in 1:ncases
        u_ret = zeros(n)
        u_ret[retained_free] = @view Ur[:, j]
        results[j] = Vector(G * u_ret)
    end
    return results
end
