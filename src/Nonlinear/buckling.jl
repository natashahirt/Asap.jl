#=
Linear Buckling Analysis
========================

Solves the eigenvalue problem: (K + λ·Kg)φ = 0

The eigenvalues λ are the buckling load factors - multipliers on the current
load level at which the structure becomes unstable.

Interpretation:
- λ > 1.0: Structure is stable under current loads
- λ < 1.0: Structure will buckle before reaching current load level
- λ = 1.0: Structure is at critical load

References:
- Cook et al. "Concepts and Applications of FEA" Ch. 14
- AISC 360-16 Appendix 7 (Direct Analysis Method)
=#

using LinearAlgebra: eigen, Symmetric

# =============================================================================
# Core Eigenvalue Solve
# =============================================================================

"""
Internal helper to solve the buckling eigenvalue problem.
Returns (load_factors, mode_shapes_reduced, n_modes).
"""
function _solve_buckling_eigen(K_ff::Matrix{Float64}, Kg_ff::Matrix{Float64}, 
                                n_requested::Int, n_dofs::Int, free_dofs::Vector{Int})
    n_free = size(K_ff, 1)
    
    # Make symmetric (handle numerical asymmetry)
    K_sym = Symmetric(0.5 * (K_ff + K_ff'))
    Kg_sym = Symmetric(0.5 * (Kg_ff + Kg_ff'))
    
    # Solve buckling eigenvalue problem: (K + λ·Kg)φ = 0
    # Reformulate as standard eigenvalue problem: K⁻¹·(-Kg)·φ = (1/λ)·φ
    # Then λ = 1/eigenvalue
    
    eigvals = Float64[]
    eigvecs = Matrix{Float64}(undef, n_free, 0)
    n = 0
    
    try
        # Form K^-1 * (-Kg)
        A = K_sym \ Matrix(-Kg_sym)
        
        # Solve standard eigenvalue problem
        eigresult = eigen(A)
        eigvals_raw = eigresult.values
        eigvecs_raw = eigresult.vectors
        
        # Convert: λ_buckling = 1/eigenvalue (for positive eigenvalues)
        valid_eigvals = Float64[]
        valid_vecs = Vector{Float64}[]
        for i in eachindex(eigvals_raw)
            μ = real(eigvals_raw[i])
            if μ > 1e-10  # Positive eigenvalue → positive buckling factor
                push!(valid_eigvals, 1.0 / μ)
                push!(valid_vecs, real.(eigvecs_raw[:, i]))
            end
        end
        
        if isempty(valid_eigvals)
            @warn "No valid buckling modes found"
            return (Float64[], zeros(n_dofs, 0), 0)
        end
        
        # Sort by buckling factor (lowest first = most critical)
        perm = sortperm(valid_eigvals)
        eigvals = valid_eigvals[perm]
        eigvecs_sorted = hcat(valid_vecs[perm]...)
        
        n = min(n_requested, length(eigvals))
        eigvals = eigvals[1:n]
        eigvecs = eigvecs_sorted[:, 1:n]
    catch e
        @warn "Buckling eigenvalue solve failed: $e"
        return (Float64[], zeros(n_dofs, 0), 0)
    end
    
    # Expand mode shapes to full DOF vector
    mode_shapes_full = zeros(n_dofs, n)
    mode_shapes_full[free_dofs, :] = eigvecs
    
    # Normalize mode shapes (unit max displacement)
    for i in 1:n
        max_val = maximum(abs.(mode_shapes_full[:, i]))
        if max_val > 1e-10
            mode_shapes_full[:, i] ./= max_val
        end
    end
    
    return (eigvals, mode_shapes_full, n)
end

# =============================================================================
# Buckling Analysis
# =============================================================================

"""
    solve_buckling!(model; n=6) -> BucklingResult

Perform linear buckling analysis to find critical load factors and mode shapes.

Solves the eigenvalue problem: (K + λ·Kg)φ = 0

# Arguments
- `model`: Structural model (FrameModel, Model, or TrussModel)
- `n::Int`: Number of buckling modes to compute (default: 6)

# Returns
- `BucklingResult`: Contains load factors, mode shapes, and summary methods

# Prerequisites
- Model must have loads applied
- Static analysis is automatically run to get internal forces for Kg

# Example
```julia
# Cantilever column under axial load
model = Model(nodes, columns, [NodeForce(top_node, [0, 0, -P])])
result = solve_buckling!(model; n=4)

# Check stability
if result.load_factors[1] > 1.0
    println("Column is stable under load P")
else
    println("Column will buckle at ", result.load_factors[1], " × P")
end

# Pretty print
print_buckling_summary(result)
```

# Notes
- Only positive eigenvalues are returned (compression buckling modes)
- Mode shapes are normalized to unit maximum
- For accuracy, use sufficient mesh refinement in slender members
"""
function solve_buckling!(model::FrameModel; n::Int = 6)
    # Ensure model is processed and solved
    if !model.processed
        process!(model)
    end
    
    # Run linear static analysis first to get internal forces
    solve!(model)
    
    # Get stiffness matrices
    K = model.S
    Kg = assemble_geometric_stiffness(model)
    
    # Extract free DOFs
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        error("No free DOFs - cannot perform buckling analysis")
    end
    
    # Limit modes to available DOFs
    n = min(n, n_free)
    
    # Extract reduced matrices (free DOFs only)
    K_ff = Matrix(K[idx, idx])
    Kg_ff = Matrix(Kg[idx, idx])
    
    eigvals, mode_shapes, n_modes = _solve_buckling_eigen(K_ff, Kg_ff, n, model.nDOFs, idx)
    
    return BucklingResult(eigvals, mode_shapes, n_modes)
end

function solve_buckling!(model::Model; n::Int = 6)
    if !model.processed
        process!(model)
    end
    
    solve!(model)
    
    K = model.S
    Kg = assemble_geometric_stiffness(model)
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        error("No free DOFs - cannot perform buckling analysis")
    end
    
    n = min(n, n_free)
    
    K_ff = Matrix(K[idx, idx])
    Kg_ff = Matrix(Kg[idx, idx])
    
    eigvals, mode_shapes, n_modes = _solve_buckling_eigen(K_ff, Kg_ff, n, model.nDOFs, idx)
    
    return BucklingResult(eigvals, mode_shapes, n_modes)
end

function solve_buckling!(model::TrussModel; n::Int = 6)
    if !model.processed
        process!(model)
    end
    
    solve!(model)
    
    K = model.S
    Kg = assemble_geometric_stiffness(model)
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        error("No free DOFs - cannot perform buckling analysis")
    end
    
    n = min(n, n_free)
    
    K_ff = Matrix(K[idx, idx])
    Kg_ff = Matrix(Kg[idx, idx])
    
    eigvals, mode_shapes, n_modes = _solve_buckling_eigen(K_ff, Kg_ff, n, model.nDOFs, idx)
    
    return BucklingResult(eigvals, mode_shapes, n_modes)
end

"""
    solve_buckling!(model::ShellModel; n=6) -> BucklingResult

Perform linear buckling analysis for shell structures.

Solves the eigenvalue problem: (K + λ·Kg)φ = 0

# Shell Buckling Physics
For plates and shells under membrane (in-plane) forces:
- Compression forces cause destabilization → positive eigenvalues
- Tension forces add stiffness → negative eigenvalues (stable)

The critical buckling load for a simply supported square plate is:
    Ncr = k·π²·D/b² where D = Et³/12(1-ν²), k=4 for uniaxial compression

# Example
```julia
# Simply supported plate under edge compression
model = ShellModel(nodes, shells, [EdgeForce(...)])
result = solve_buckling!(model; n=4)

if result.load_factors[1] < 1.0
    println("Plate will buckle before reaching applied load")
end
```
"""
function solve_buckling!(model::ShellModel; n::Int = 6)
    if !model.processed
        process!(model)
    end
    
    solve!(model)
    
    K = model.S
    Kg = assemble_geometric_stiffness(model)
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        error("No free DOFs - cannot perform buckling analysis")
    end
    
    n = min(n, n_free)
    
    K_ff = Matrix(K[idx, idx])
    Kg_ff = Matrix(Kg[idx, idx])
    
    eigvals, mode_shapes, n_modes = _solve_buckling_eigen(K_ff, Kg_ff, n, model.nDOFs, idx)
    
    return BucklingResult(eigvals, mode_shapes, n_modes)
end

"""
    solve_buckling!(model::ShellModel, σ_uniform::Vector{Float64}; n=6) -> BucklingResult

Perform linear buckling analysis with prescribed uniform membrane forces.

This is useful for validation against analytical solutions where you know
the stress state (e.g., uniform compression).

# Arguments
- `model::ShellModel`: Processed shell model (must have BCs applied via node.dof)
- `σ_uniform::Vector{Float64}`: [Nxx, Nyy, Nxy] membrane forces [N/m]
  - Nxx < 0 for compression in x
  - Nyy < 0 for compression in y
- `n::Int`: Number of buckling modes to compute (default: 6)

# Returns
- `BucklingResult`: Contains load factors, mode shapes, and summary methods

# Example
```julia
# Simply supported square plate under uniform uniaxial compression
# Analytical Ncr = 4π²D/a² for square plate

D = E * t^3 / (12 * (1 - ν^2))
Ncr_analytical = 4 * π^2 * D / a^2

# Apply unit compression and find load factor
σ_unit = [-1.0, 0.0, 0.0]  # Unit compression in x
result = solve_buckling!(model, σ_unit; n=4)

# Critical load is λ * |σ_unit|
Ncr_computed = result.load_factors[1] * 1.0
```
"""
function solve_buckling!(model::ShellModel, σ_uniform::Vector{Float64}; n::Int = 6)
    if !model.processed
        process!(model)
    end
    
    K = model.S
    Kg = assemble_geometric_stiffness(model, σ_uniform)
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        error("No free DOFs - cannot perform buckling analysis")
    end
    
    n = min(n, n_free)
    
    K_ff = Matrix(K[idx, idx])
    Kg_ff = Matrix(Kg[idx, idx])
    
    eigvals, mode_shapes, n_modes = _solve_buckling_eigen(K_ff, Kg_ff, n, model.nDOFs, idx)
    
    return BucklingResult(eigvals, mode_shapes, n_modes)
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    critical_load_factor(model) -> Float64

Convenience function to get just the first (most critical) buckling load factor.
"""
function critical_load_factor(model::Union{FrameModel, Model, TrussModel, ShellModel})
    result = solve_buckling!(model; n=1)
    return result.n_modes > 0 ? result.load_factors[1] : Inf
end

"""
    is_stable(model) -> Bool

Check if structure is stable under current loading (λ_cr > 1.0).
"""
function is_stable(model::Union{FrameModel, Model, TrussModel, ShellModel})
    return critical_load_factor(model) > 1.0
end
