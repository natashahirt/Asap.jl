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

