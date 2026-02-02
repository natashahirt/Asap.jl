#=
Modal Analysis Module for Asap.jl
=================================

Provides eigenvalue analysis for natural frequencies and mode shapes.

Mathematical basis:
  (K - ω²M)φ = 0

where K is stiffness matrix, M is mass matrix, ω is natural frequency,
and φ is the mode shape vector.

References:
- Cook, Malkus, Plesha "Concepts and Applications of FEA"
- FinEtoolsFlexStructures.jl by Petr Krysl - MIT License
=#

using LinearAlgebra: eigen, Symmetric, I
using SparseArrays: sparse, SparseMatrixCSC

export modal_analysis, natural_frequencies, mode_shapes
export assemble_mass_matrix, assemble_mass_matrix!

# =============================================================================
# Mass Matrix Assembly - FrameModel
# =============================================================================

"""
    assemble_mass_matrix(model; type=MASS_CONSISTENT)

Assemble the global mass matrix for a structural model.

# Arguments
- `model`: Processed structural model (FrameModel, ShellModel, or Model)
- `type::MassMatrixType`: Mass matrix formulation (default: `MASS_CONSISTENT`)

# Returns
- `SparseMatrixCSC{Float64}`: Global mass matrix
"""
function assemble_mass_matrix(model::FrameModel; type::MassMatrixType = MASS_CONSISTENT)
    n_dof = model.nDOFs
    M = spzeros(Float64, n_dof, n_dof)
    assemble_mass_matrix!(M, model; type=type)
    return M
end

function assemble_mass_matrix!(M::SparseMatrixCSC{Float64}, model::FrameModel; 
                               type::MassMatrixType = MASS_CONSISTENT)
    for element in model.elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    return M
end

function assemble_mass_matrix!(M::Matrix{Float64}, model::FrameModel; 
                               type::MassMatrixType = MASS_CONSISTENT)
    for element in model.elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    return M
end

# =============================================================================
# Mass Matrix Assembly - ShellModel
# =============================================================================

function assemble_mass_matrix(model::ShellModel; type::MassMatrixType = MASS_CONSISTENT)
    n_dof = model.nDOFs
    M = spzeros(Float64, n_dof, n_dof)
    assemble_mass_matrix!(M, model; type=type)
    return M
end

function assemble_mass_matrix!(M::SparseMatrixCSC{Float64}, model::ShellModel; 
                               type::MassMatrixType = MASS_CONSISTENT)
    for element in model.elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    return M
end

function assemble_mass_matrix!(M::Matrix{Float64}, model::ShellModel; 
                               type::MassMatrixType = MASS_CONSISTENT)
    for element in model.elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    return M
end

# =============================================================================
# Mass Matrix Assembly - Unified Model
# =============================================================================

function assemble_mass_matrix(model::Model; type::MassMatrixType = MASS_CONSISTENT)
    n_dof = model.nDOFs
    M = spzeros(Float64, n_dof, n_dof)
    assemble_mass_matrix!(M, model; type=type)
    return M
end

function assemble_mass_matrix!(M::SparseMatrixCSC{Float64}, model::Model; 
                               type::MassMatrixType = MASS_CONSISTENT)
    # Frame elements
    for element in model.frame_elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    
    # Shell elements
    for element in model.shell_elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    
    return M
end

function assemble_mass_matrix!(M::Matrix{Float64}, model::Model; 
                               type::MassMatrixType = MASS_CONSISTENT)
    # Frame elements
    for element in model.frame_elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    
    # Shell elements
    for element in model.shell_elements
        M_global = global_mass(element; type=type)
        gid = element.globalID
        
        for (i, gi) in enumerate(gid)
            for (j, gj) in enumerate(gid)
                M[gi, gj] += M_global[i, j]
            end
        end
    end
    
    return M
end

# =============================================================================
# Modal Analysis Result
# =============================================================================

"""
    ModalResult

Result container for modal analysis.

# Fields
- `frequencies::Vector{Float64}`: Natural frequencies [Hz]
- `omegas::Vector{Float64}`: Angular frequencies [rad/s]
- `periods::Vector{Float64}`: Periods [s]
- `mode_shapes::Matrix{Float64}`: Mode shape vectors (columns)
- `n_modes::Int`: Number of computed modes
- `mass_type::MassMatrixType`: Mass formulation used
"""
struct ModalResult
    frequencies::Vector{Float64}   # Natural frequencies [Hz]
    omegas::Vector{Float64}        # Angular frequencies [rad/s]
    periods::Vector{Float64}       # Periods [s]
    mode_shapes::Matrix{Float64}   # Mode shapes (each column is a mode)
    n_modes::Int                   # Number of modes computed
    mass_type::MassMatrixType      # Mass formulation used
end

# =============================================================================
# Modal Analysis - Generic Implementation
# =============================================================================

"""
    _modal_analysis_impl(S, nDOFs, freeDOFs, model, n_modes, mass_type)

Internal implementation of modal analysis.
"""
function _modal_analysis_impl(S::SparseMatrixCSC, nDOFs::Int, freeDOFs::Vector{Int},
                              model, n_modes::Int, mass_type::MassMatrixType)
    
    n_free = length(freeDOFs)
    
    if n_free == 0
        error("No free DOFs in model - cannot perform modal analysis")
    end
    
    # Limit modes to available DOFs
    n_modes = min(n_modes, n_free)
    
    # Extract reduced stiffness matrix (free DOFs only)
    K_full = Matrix(S)
    K_ff = K_full[freeDOFs, freeDOFs]
    
    # Assemble and reduce mass matrix
    M_full = zeros(nDOFs, nDOFs)
    assemble_mass_matrix!(M_full, model; type=mass_type)
    M_ff = M_full[freeDOFs, freeDOFs]
    
    # Solve generalized eigenvalue problem: K*φ = ω²*M*φ
    # Make matrices symmetric (numerical symmetry)
    K_sym = Symmetric(0.5 * (K_ff + K_ff'))
    M_sym = Symmetric(0.5 * (M_ff + M_ff'))
    
    # Solve generalized eigenvalue problem
    eigvals, eigvecs = eigen(K_sym, M_sym)
    
    # Filter out negative/zero eigenvalues (rigid body modes or numerical noise)
    valid_idx = findall(λ -> λ > 1e-10, eigvals)
    
    if length(valid_idx) < n_modes
        @warn "Only $(length(valid_idx)) valid modes found (requested $n_modes)"
        n_modes = length(valid_idx)
    end
    
    # Take first n_modes valid eigenvalues
    mode_idx = valid_idx[1:n_modes]
    
    omega_sq = eigvals[mode_idx]
    omegas = sqrt.(omega_sq)               # Angular frequencies [rad/s]
    frequencies = omegas ./ (2π)           # Natural frequencies [Hz]
    periods = 1.0 ./ frequencies           # Periods [s]
    
    # Extract mode shapes and expand to full DOF vector
    mode_shapes_reduced = eigvecs[:, mode_idx]
    mode_shapes_full = zeros(nDOFs, n_modes)
    mode_shapes_full[freeDOFs, :] = mode_shapes_reduced
    
    # Normalize mode shapes (mass normalization: φᵀMφ = 1)
    for i in 1:n_modes
        φ = mode_shapes_full[:, i]
        mass_norm = sqrt(φ' * M_full * φ)
        if mass_norm > 1e-10
            mode_shapes_full[:, i] ./= mass_norm
        end
    end
    
    return ModalResult(frequencies, omegas, periods, mode_shapes_full, n_modes, mass_type)
end

# =============================================================================
# Modal Analysis - Public API
# =============================================================================

"""
    modal_analysis(model; n_modes=6, mass_type=MASS_CONSISTENT)

Perform eigenvalue analysis to find natural frequencies and mode shapes.

Works with FrameModel, ShellModel, or unified Model.

# Arguments
- `model`: Processed structural model
- `n_modes::Int`: Number of modes to compute (default: 6)
- `mass_type::MassMatrixType`: Mass matrix formulation

# Returns
- `ModalResult`: Container with frequencies, mode shapes, etc.

# Example
```julia
process!(model)
result = modal_analysis(model; n_modes=10)
println("First natural frequency: \$(result.frequencies[1]) Hz")
```
"""
function modal_analysis(model::FrameModel; n_modes::Int = 6, 
                        mass_type::MassMatrixType = MASS_CONSISTENT)
    if !model.processed
        error("Model must be processed before modal analysis. Call process!(model) first.")
    end
    return _modal_analysis_impl(model.S, model.nDOFs, model.freeDOFs, model, n_modes, mass_type)
end

function modal_analysis(model::ShellModel; n_modes::Int = 6, 
                        mass_type::MassMatrixType = MASS_CONSISTENT)
    if !model.processed
        error("Model must be processed before modal analysis. Call process!(model) first.")
    end
    return _modal_analysis_impl(model.S, model.nDOFs, model.freeDOFs, model, n_modes, mass_type)
end

function modal_analysis(model::Model; n_modes::Int = 6, 
                        mass_type::MassMatrixType = MASS_CONSISTENT)
    if !model.processed
        error("Model must be processed before modal analysis. Call process!(model) first.")
    end
    return _modal_analysis_impl(model.S, model.nDOFs, model.freeDOFs, model, n_modes, mass_type)
end

# =============================================================================
# Convenience Functions
# =============================================================================

"""
    natural_frequencies(model; n_modes=6, mass_type=MASS_CONSISTENT)

Convenience function to get just the natural frequencies [Hz].
"""
function natural_frequencies(model::Union{FrameModel, ShellModel, Model}; 
                            n_modes::Int = 6,
                            mass_type::MassMatrixType = MASS_CONSISTENT)
    result = modal_analysis(model; n_modes=n_modes, mass_type=mass_type)
    return result.frequencies
end

"""
    mode_shapes(model; n_modes=6, mass_type=MASS_CONSISTENT)

Convenience function to get just the mode shape vectors.
"""
function mode_shapes(model::Union{FrameModel, ShellModel, Model}; 
                    n_modes::Int = 6,
                    mass_type::MassMatrixType = MASS_CONSISTENT)
    result = modal_analysis(model; n_modes=n_modes, mass_type=mass_type)
    return result.mode_shapes
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    print_modal_summary(result::ModalResult)

Print a formatted summary of modal analysis results.
"""
function print_modal_summary(result::ModalResult)
    println("Modal Analysis Summary")
    println("=" ^ 50)
    println("Number of modes: $(result.n_modes)")
    println("Mass formulation: $(result.mass_type)")
    println()
    println("Mode |  Frequency (Hz)  |  Period (s)  | ω (rad/s)")
    println("-" ^ 50)
    for i in 1:result.n_modes
        f_str = lpad(round(result.frequencies[i], digits=4), 14)
        p_str = lpad(round(result.periods[i], digits=4), 10)
        w_str = lpad(round(result.omegas[i], digits=4), 10)
        println(" $(lpad(i, 3)) | $f_str   | $p_str   | $w_str")
    end
end

# =============================================================================
# Primary Modal API: modal! and modal
# =============================================================================

"""
    modal!(model; n=6, mass_type=MASS_CONSISTENT)

Perform modal analysis on a structural model (mutating).

Automatically processes the model if not already processed.
This is the recommended way to run modal analysis.

# Arguments
- `model`: Structural model (Model, FrameModel, ShellModel)
- `n::Int`: Number of modes to compute (default: 6)
- `mass_type::MassMatrixType`: Mass formulation (default: MASS_CONSISTENT)

# Returns
- `ModalResult`: Container with frequencies, periods, mode shapes

# Example
```julia
model = Model(nodes, elements, loads)
result = modal!(model; n=10)

# Access results
result.frequencies   # Natural frequencies [Hz]
result.periods       # Periods [s]
result.mode_shapes   # Mode shape matrix (columns are modes)

# Pretty print
print_modal_summary(result)
```
"""
function modal!(model::FrameModel; n::Int = 6, mass_type::MassMatrixType = MASS_CONSISTENT)
    if !model.processed
        process!(model)
    end
    return modal_analysis(model; n_modes=n, mass_type=mass_type)
end

function modal!(model::ShellModel; n::Int = 6, mass_type::MassMatrixType = MASS_CONSISTENT)
    if !model.processed
        process!(model)
    end
    return modal_analysis(model; n_modes=n, mass_type=mass_type)
end

function modal!(model::Model; n::Int = 6, mass_type::MassMatrixType = MASS_CONSISTENT)
    if !model.processed
        process!(model)
    end
    return modal_analysis(model; n_modes=n, mass_type=mass_type)
end

"""
    modal(model; n=6, mass_type=MASS_CONSISTENT)

Perform modal analysis on a structural model (non-mutating).

Returns the modal result without modifying the model. Model must already
be processed.

# Arguments
- `model`: Processed structural model
- `n::Int`: Number of modes to compute (default: 6)
- `mass_type::MassMatrixType`: Mass formulation (default: MASS_CONSISTENT)

# Returns
- `ModalResult`: Container with frequencies, periods, mode shapes

# Example
```julia
solve!(model)  # Process and solve statics first
result = modal(model; n=6)
```
"""
function modal(model::Union{FrameModel, ShellModel, Model}; n::Int = 6, 
               mass_type::MassMatrixType = MASS_CONSISTENT)
    return modal_analysis(model; n_modes=n, mass_type=mass_type)
end