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

