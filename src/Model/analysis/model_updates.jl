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
  sparse reconstruction), and the next solve refactors `S` reusing the
  previous symbolic analysis.  Use when only E, t, ν, I, A, or load
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
        _invalidate_factor_values!(model)
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
        _recompute_shell_K!(model.shell_elements)
        populate_loads!(model)
        _update_S_values!(model)
        _invalidate_factor_values!(model)
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
        _recompute_shell_K!(model.elements)
        populate_loads!(model)
        _update_S_values!(model)
        _invalidate_factor_values!(model)
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
        _invalidate_factor_values!(model)
    else
        process_elements!(model)
        populate_loads!(model)
        create_S!(model)
        model._factorization = nothing
    end

    model.processed = true
end


