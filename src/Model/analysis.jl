# =============================================================================
# Analysis Types
# =============================================================================

"""
Supported analysis types for `solve!(model, analysis_type)`.

**Single Analyses:**
- `:static` - Linear static (K·u = F)
- `:modal` - Modal/eigenvalue analysis
- `:buckling` - Linear buckling eigenvalue
- `:pdelta` - P-delta iteration with geometric stiffness
- `:nonlinear` - Full Newton-Raphson nonlinear static

**Composite Analyses:**
- `:stability` - Buckling + P-delta (returns NamedTuple)
- `:all` - All valid analyses for model type (returns NamedTuple)
"""
const ANALYSIS_TYPES = (
    :static, :modal, :buckling, :pdelta, :nonlinear,  # Single
    :stability, :all                                    # Composite
)

"""
    available_analyses(model) -> Tuple{Symbol...}

Return the analysis types that are valid for this model type.

- `FrameModel`: All analyses (static, modal, buckling, pdelta, nonlinear)
- `ShellModel`: Static and modal only (no geometric stiffness for shells alone)
- `TrussModel`: Static and buckling only (no mass matrix for modal)
- `Model`: Depends on elements - full if has frames, limited if shells only
"""
available_analyses(::FrameModel) = (:static, :modal, :buckling, :pdelta, :nonlinear)
available_analyses(::ShellModel) = (:static, :modal)
available_analyses(::TrussModel) = (:static, :buckling)

function available_analyses(model::Model)
    if !isempty(model.frame_elements)
        # Has frame elements - full analysis suite
        return (:static, :modal, :buckling, :pdelta, :nonlinear)
    else
        # Shell-only model
        return (:static, :modal)
    end
end

"""
    supports_analysis(model, analysis::Symbol) -> Bool

Check if a model supports a specific analysis type.
"""
function supports_analysis(model::AbstractModel, analysis::Symbol)
    analysis ∈ (:stability, :all) && return true  # Composite handled separately
    return analysis ∈ available_analyses(model)
end

# =============================================================================
# Model Processing
# =============================================================================

"""
    process!(model::FrameModel)

Process a structural model: add linkages between nodes and elements, 
determine DOF orders, generate the load vectors P, Pf, and assemble 
the global stiffness matrix, S.
"""
function process!(model::FrameModel)
    make_ids!(model)

    if any(e -> e isa BridgeElement, model.elements)
        processBridge!(model)
        make_ids!(model)
    else
        process_elements!(model)
    end

    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model._factorization = nothing
    model._elemental_loads = nothing
    model.processed = true
end

"""
    process!(model::ShellModel)

Process a shell model: add linkages between nodes and elements,
determine DOF orders, generate the load vector P, and assemble
the global stiffness matrix, S.
"""
function process!(model::ShellModel)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model._factorization = nothing
    model.processed = true
end

"""
    process!(model::Model)

Process a unified model: add linkages between nodes and elements,
determine DOF orders, generate the load vectors P, Pf, and assemble
the global stiffness matrix, S.

Works with mixed frame+shell models or single-type models.
"""
function process!(model::Model)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model._factorization = nothing
    model._elemental_loads = nothing
    model.processed = true
end

"""
    process!(model::TrussModel)

Process a structural truss model.
"""
function process!(model::TrussModel)
    make_ids!(model)
    process_elements!(model)
    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

    model._factorization = nothing
    model.processed = true
end

# =============================================================================
# Solving
# =============================================================================

"""
    solve!(model::FrameModel; reprocess = false)

Solve for the nodal displacements of a structural model. 
`reprocess = true` reevaluates all node/element properties and reassembles the global stiffness matrix.
"""
function solve!(model::FrameModel; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        for element in model.elements
            if element isa Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    fact = _get_factorization(model)
    U = fact \ F

    model.compliance = dot(U, F)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

"""
    solve!(model::ShellModel; reprocess = false)

Solve for the nodal displacements of a shell model.
"""
function solve!(model::ShellModel; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    P_free = @view model.P[idx]
    fact = _get_factorization(model)
    U = fact \ P_free

    model.compliance = dot(U, P_free)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

"""
    solve!(model::Model; reprocess = false)

Solve for the nodal displacements of a unified model.
Works with mixed frame+shell models or single-type models.
"""
function solve!(model::Model; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        for element in model.frame_elements
            if element isa Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    fact = _get_factorization(model)
    U = fact \ F

    model.compliance = dot(U, F)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

"""
    solve!(model::TrussModel; reprocess = false)

Solve for the nodal displacements of a structural truss model.
"""
function solve!(model::TrussModel; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    P_free = @view model.P[idx]
    fact = _get_factorization(model)
    U = fact \ P_free

    model.compliance = dot(U, P_free)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

# =============================================================================
# Solving with Custom Loads
# =============================================================================

"""
    solve(model::FrameModel, L::Vector{AbstractLoad})

Return the displacement vector under a given load set `L`, reusing the cached
stiffness factorization.
"""
function solve(model::FrameModel, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    fact = _get_factorization(model)
    
    u = zeros(model.nDOFs)
    u[idx] = fact \ F[idx]
    return u
end

"""
    solve(model::Model, L::Vector{AbstractLoad})

Return the displacement vector under a given load set `L`, reusing the cached
stiffness factorization.
"""
function solve(model::Model, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    fact = _get_factorization(model)
    
    u = zeros(model.nDOFs)
    u[idx] = fact \ F[idx]
    return u
end

"""
    solve(model::FrameModel, cases::Vector{<:Vector{<:AbstractLoad}})

Solve multiple load cases against the same stiffness matrix. Returns a vector
of displacement vectors, one per case. Uses a single factorization and a
multi-RHS LAPACK solve for efficiency.
"""
function solve(model::FrameModel, cases::Vector{<:Vector{<:AbstractLoad}})
    model.processed || process!(model)
    fact = _get_factorization(model)
    idx = model.freeDOFs
    n = model.nDOFs
    ncases = length(cases)

    # Build force matrix — each column is a load case
    F = Matrix{Float64}(undef, length(idx), ncases)
    for (j, loads) in enumerate(cases)
        Fj = create_F(model, loads)
        F[:, j] = Fj[idx]
    end

    # Multi-RHS solve (LAPACK level-3 BLAS)
    U = fact \ F

    # Unpack into full DOF vectors
    results = Vector{Vector{Float64}}(undef, ncases)
    for j in 1:ncases
        u = zeros(n)
        u[idx] = @view U[:, j]
        results[j] = u
    end
    return results
end

"""
    solve(model::Model, cases::Vector{<:Vector{<:AbstractLoad}})

Solve multiple load cases against the same stiffness matrix. Returns a vector
of displacement vectors, one per case. Uses a single factorization and a
multi-RHS LAPACK solve for efficiency.
"""
function solve(model::Model, cases::Vector{<:Vector{<:AbstractLoad}})
    model.processed || process!(model)
    fact = _get_factorization(model)
    idx = model.freeDOFs
    n = model.nDOFs
    ncases = length(cases)

    # Build force matrix — each column is a load case
    F = Matrix{Float64}(undef, length(idx), ncases)
    for (j, loads) in enumerate(cases)
        Fj = create_F(model, loads)
        F[:, j] = Fj[idx]
    end

    # Multi-RHS solve (LAPACK level-3 BLAS)
    U = fact \ F

    # Unpack into full DOF vectors
    results = Vector{Vector{Float64}}(undef, ncases)
    for j in 1:ncases
        u = zeros(n)
        u[idx] = @view U[:, j]
        results[j] = u
    end
    return results
end

"""
    solve!(model::FrameModel, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::FrameModel, L::Vector{<:AbstractLoad})
    model.loads = L
    model._elemental_loads = nothing
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve!(model::Model, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::Model, L::Vector{<:AbstractLoad})
    model.loads = L
    model._elemental_loads = nothing
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector under a given load set `L`, reusing the cached
stiffness factorization.
"""
function solve(model::TrussModel, L::Vector{NodeForce})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    fact = _get_factorization(model)
    
    u = zeros(model.nDOFs)
    u[idx] = fact \ F[idx]
    return u
end

"""
    solve!(model::TrussModel, L::Vector{NodeForce})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::TrussModel, L::Vector{NodeForce})
    model.loads = L
    process!(model)
    solve!(model)
    post_process!(model)
end

# =============================================================================
# Unified Analysis API
# =============================================================================

"""
    solve!(model, analysis::Symbol; kwargs...) -> Result

Unified analysis entry point. Run a specific analysis type on the model.

# Analysis Types

**Single Analyses:**
- `:static` - Linear static analysis (default, K·u = F)
- `:modal` - Modal analysis (natural frequencies and mode shapes)
- `:buckling` - Linear buckling (critical load factors)
- `:pdelta` - P-delta iteration (second-order effects)
- `:nonlinear` - Nonlinear static (Newton-Raphson)

**Composite Analyses:**
- `:stability` - Buckling + P-delta together
- `:all` - All analyses

# Examples
```julia
# Linear static (returns nothing, updates model)
solve!(model, :static)

# Modal analysis
result = solve!(model, :modal; n=6)
println(result.frequencies)

# Buckling
result = solve!(model, :buckling; n=4)
println("Critical λ = ", result.load_factors[1])

# P-delta
result = solve!(model, :pdelta; max_iter=10)
println("Amplification = ", result.amplification)

# Stability check (buckling + pdelta)
results = solve!(model, :stability)
println(results.buckling.load_factors[1])
println(results.pdelta.amplification)

# Everything
results = solve!(model, :all)
```

# Keyword Arguments
- `:modal`: `n` (number of modes), `mass_type`
- `:buckling`: `n` (number of modes)
- `:pdelta`: `max_iter`, `tol`, `verbose`
- `:nonlinear`: `n_steps`, `max_iter`, `tol`, `verbose`
"""
function solve!(model::AbstractModel, analysis::Symbol; kwargs...)
    analysis ∈ ANALYSIS_TYPES || error("Unknown analysis type: $analysis. Valid types: $ANALYSIS_TYPES")
    
    # Check if analysis is supported for this model type (except composites)
    if analysis ∉ (:stability, :all) && !supports_analysis(model, analysis)
        avail = available_analyses(model)
        error("Analysis :$analysis not supported for $(typeof(model)). Available: $avail")
    end
    
    # Single analyses
    if analysis == :static
        _solve_static!(model; kwargs...)
        return nothing
        
    elseif analysis == :modal
        return _solve_modal!(model; kwargs...)
        
    elseif analysis == :buckling
        return _solve_buckling!(model; kwargs...)
        
    elseif analysis == :pdelta
        return _solve_pdelta!(model; kwargs...)
        
    elseif analysis == :nonlinear
        return _solve_nonlinear!(model; kwargs...)
        
    # Composite analyses
    elseif analysis == :stability
        # Stability requires buckling + pdelta (frame elements needed)
        avail = available_analyses(model)
        if :buckling ∉ avail || :pdelta ∉ avail
            error(":stability requires frame elements. Available for this model: $avail")
        end
        n_buck = get(kwargs, :n, 4)
        max_iter = get(kwargs, :max_iter, 10)
        tol = get(kwargs, :tol, 1e-3)
        buck = _solve_buckling!(model; n=n_buck)
        pdelta = _solve_pdelta!(model; max_iter=max_iter, tol=tol)
        return (buckling=buck, pdelta=pdelta)
        
    elseif analysis == :all
        # Run all analyses that are valid for this model type
        avail = available_analyses(model)
        n_modal = get(kwargs, :n_modal, 6)
        n_buck = get(kwargs, :n_buckling, 4)
        
        # Always run static first
        _solve_static!(model)
        
        # Build results based on what's available
        results = Dict{Symbol, Any}()
        
        if :modal ∈ avail
            results[:modal] = _solve_modal!(model; n=n_modal)
        end
        if :buckling ∈ avail
            results[:buckling] = _solve_buckling!(model; n=n_buck)
        end
        if :pdelta ∈ avail
            results[:pdelta] = _solve_pdelta!(model)
        end
        
        # Convert to NamedTuple for nice access
        return NamedTuple(results)
    end
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
function _get_factorization(model::AbstractModel)
    if model._factorization !== nothing
        return model._factorization
    end
    idx = model.freeDOFs
    K = Symmetric(model.S[idx, idx])
    fact = cholesky(K; check=false)
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

# Static - wrapper for existing solve!
function _solve_static!(model::FrameModel; reprocess=false, postprocess::Symbol=:all)
    if !model.processed || reprocess
        for element in model.elements
            if element isa Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    fact = _get_factorization(model)
    U = fact \ F

    model.compliance = dot(U, F)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

function _solve_static!(model::ShellModel; reprocess=false, postprocess::Symbol=:all)
    if !model.processed || reprocess
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    P_free = @view model.P[idx]
    fact = _get_factorization(model)
    U = fact \ P_free

    model.compliance = dot(U, P_free)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

function _solve_static!(model::Model; reprocess=false, postprocess::Symbol=:all)
    if !model.processed || reprocess
        for element in model.frame_elements
            if element isa Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    fact = _get_factorization(model)
    U = fact \ F

    model.compliance = dot(U, F)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

function _solve_static!(model::TrussModel; reprocess=false, postprocess::Symbol=:all)
    if !model.processed || reprocess
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    P_free = @view model.P[idx]
    fact = _get_factorization(model)
    U = fact \ P_free

    model.compliance = dot(U, P_free)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

# Modal - delegates to modal!
function _solve_modal!(model; n::Int=6, mass_type::MassMatrixType=MASS_CONSISTENT)
    return modal!(model; n=n, mass_type=mass_type)
end

# Buckling - delegates to solve_buckling!
function _solve_buckling!(model; n::Int=6)
    return solve_buckling!(model; n=n)
end

# P-delta - delegates to solve_pdelta!
function _solve_pdelta!(model; max_iter::Int=10, tol::Float64=1e-3, verbose::Bool=false)
    return solve_pdelta!(model; max_iter=max_iter, tol=tol, verbose=verbose)
end

# Nonlinear - delegates to solve_nonlinear!
function _solve_nonlinear!(model; n_steps::Int=10, max_iter::Int=20, tol::Float64=1e-4, verbose::Bool=false)
    return solve_nonlinear!(model; n_steps=n_steps, max_iter=max_iter, tol=tol, verbose=verbose)
end

# =============================================================================
# update!  —  public API for property/load changes (topology unchanged)
# =============================================================================

"""
    update!(model; values_only=false)

Rebuild element stiffness matrices, load vector, and global stiffness
matrix after property or load changes.  Topology must be unchanged
(same nodes, elements, DOFs) — use `process!` after topology changes.

# Keyword arguments
- `values_only=false` (default): full property update — recomputes element
  geometry (LCS, R, length), stiffness K, loads P, and reassembles S from
  COO triplets.  Use after node positions moved or section connectivity
  changed.

- `values_only=true`: lightweight update — only recomputes element K
  matrices from current section/material properties and rebuilds P.
  Updates S values **in-place** (same sparsity pattern, avoids COO →
  sparse reconstruction).  Use when only E, t, ν, I, A, or load
  magnitudes changed but geometry is unchanged.
"""
function update!(model::FrameModel; values_only::Bool = false, loads_only::Bool = false)
    for element in model.elements
        element isa Element && fill!(element.Q, 0.0)
    end

    if loads_only
        # Only load magnitudes changed — K and S untouched, keep factorization
        populate_loads!(model)
    elseif values_only
        # K-only pass: skip lcs!, R!, length! — geometry unchanged
        for element in model.elements
            global_K!(element)
        end
        populate_loads!(model)
        _update_S_values!(model)
        model._factorization = nothing
    else
        # Full element reprocess + S rebuild
        process_elements!(model)
        populate_loads!(model)
        create_S!(model)
        model._factorization = nothing
    end

    model._elemental_loads = nothing
    model.processed = true
end

function update!(model::Model; values_only::Bool = false, loads_only::Bool = false)
    for element in model.frame_elements
        element isa Element && fill!(element.Q, 0.0)
    end

    if loads_only
        # Only load magnitudes changed — K and S untouched, keep factorization
        populate_loads!(model)
    elseif values_only
        # K-only pass: skip geometry — reuse existing LCS, R, length, area
        for element in model.frame_elements
            global_K!(element)
        end
        # Shell global_K! is thread-safe (writes only to per-element fields)
        shells = model.shell_elements
        Threads.@threads for i in eachindex(shells)
            global_K!(shells[i])
        end
        populate_loads!(model)
        _update_S_values!(model)
        model._factorization = nothing
    else
        # Full element reprocess + S rebuild
        process_elements!(model)
        populate_loads!(model)
        create_S!(model)
        model._factorization = nothing
    end

    model._elemental_loads = nothing
    model.processed = true
end

function update!(model::ShellModel; values_only::Bool = false, loads_only::Bool = false)
    if loads_only
        # Only load magnitudes changed — K and S untouched, keep factorization
        populate_loads!(model)
    elseif values_only
        # Shell global_K! is thread-safe (writes only to per-element fields)
        elems = model.elements
        Threads.@threads for i in eachindex(elems)
            global_K!(elems[i])
        end
        populate_loads!(model)
        _update_S_values!(model)
        model._factorization = nothing
    else
        process_elements!(model)
        populate_loads!(model)
        create_S!(model)
        model._factorization = nothing
    end

    model.processed = true
end

function update!(model::TrussModel; values_only::Bool = false, loads_only::Bool = false)
    if loads_only
        # Only load magnitudes changed — K and S untouched, keep factorization
        populate_loads!(model)
    elseif values_only
        for element in model.elements
            global_K!(element)
        end
        populate_loads!(model)
        _update_S_values!(model)
        model._factorization = nothing
    else
        process_elements!(model)
        populate_loads!(model)
        create_S!(model)
        model._factorization = nothing
    end

    model.processed = true
end


# =============================================================================
# _update_S_values!  —  in-place stiffness scatter (same sparsity pattern)
# =============================================================================

"""
    _update_S_values!(model)

Zero the existing sparse matrix `model.S` and scatter updated element K
values back into it.  Reuses the CSC sparsity pattern from the initial
`create_S!` call — avoids COO allocation + `sparse()` reconstruction.
"""
function _update_S_values!(model::FrameModel)
    fill!(model.S.nzval, 0.0)
    for element in model.elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i in 1:n, j in 1:n
            model.S[idx[i], idx[j]] += element.K[i, j]
        end
    end
end

function _update_S_values!(model::ShellModel)
    fill!(model.S.nzval, 0.0)
    for element in model.elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i in 1:n, j in 1:n
            model.S[idx[i], idx[j]] += element.K[i, j]
        end
    end
end

function _update_S_values!(model::TrussModel)
    fill!(model.S.nzval, 0.0)
    for element in model.elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i in 1:n, j in 1:n
            model.S[idx[i], idx[j]] += element.K[i, j]
        end
    end
end

function _update_S_values!(model::Model)
    fill!(model.S.nzval, 0.0)
    for element in model.frame_elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i in 1:n, j in 1:n
            model.S[idx[i], idx[j]] += element.K[i, j]
        end
    end
    for element in model.shell_elements
        idx = element.globalID
        n = length(idx)
        @inbounds for i in 1:n, j in 1:n
            model.S[idx[i], idx[j]] += element.K[i, j]
        end
    end
end
