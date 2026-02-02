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

    if any(typeof.(model.elements) .<: BridgeElement)
        processBridge!(model)
        make_ids!(model)
    else
        process_elements!(model)
    end

    populate_DOF_indices!(model)
    populate_loads!(model)
    create_S!(model)

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
function solve!(model::FrameModel; reprocess = false)
    if !model.processed || reprocess
        for element in model.elements
            if typeof(element) <: Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    model.compliance = U' * F
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

"""
    solve!(model::ShellModel; reprocess = false)

Solve for the nodal displacements of a shell model.
"""
function solve!(model::ShellModel; reprocess = false)
    if !model.processed || reprocess
        process!(model)
    end

    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    model.compliance = U' * model.P[idx]
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

"""
    solve!(model::Model; reprocess = false)

Solve for the nodal displacements of a unified model.
Works with mixed frame+shell models or single-type models.
"""
function solve!(model::Model; reprocess = false)
    if !model.processed || reprocess
        # Clear fixed-end forces for frame elements
        for element in model.frame_elements
            if typeof(element) <: Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    model.compliance = U' * F
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

"""
    solve!(model::TrussModel; reprocess = false)

Solve for the nodal displacements of a structural truss model.
"""
function solve!(model::TrussModel; reprocess = false)
    if !model.processed || reprocess
        process!(model)
    end

    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    model.compliance = U' * model.P[idx]
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

# =============================================================================
# Solving with Custom Loads
# =============================================================================

"""
    solve(model::FrameModel, L::Vector{AbstractLoad})

Return the displacement vector under a given load set L.
"""
function solve(model::FrameModel, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve(model::Model, L::Vector{AbstractLoad})

Return the displacement vector under a given load set L.
"""
function solve(model::Model, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::FrameModel, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::FrameModel, L::Vector{<:AbstractLoad})
    model.loads = L
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
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector under a given load set L.
"""
function solve(model::TrussModel, L::Vector{NodeForce})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

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

# Static - wrapper for existing solve!
function _solve_static!(model::FrameModel; reprocess=false)
    if !model.processed || reprocess
        for element in model.elements
            if typeof(element) <: Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    model.compliance = U' * F
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

function _solve_static!(model::ShellModel; reprocess=false)
    if !model.processed || reprocess
        process!(model)
    end

    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    model.compliance = U' * model.P[idx]
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

function _solve_static!(model::Model; reprocess=false)
    if !model.processed || reprocess
        for element in model.frame_elements
            if typeof(element) <: Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    model.compliance = U' * F
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
end

function _solve_static!(model::TrussModel; reprocess=false)
    if !model.processed || reprocess
        process!(model)
    end

    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    model.compliance = U' * model.P[idx]
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    post_process!(model)
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
