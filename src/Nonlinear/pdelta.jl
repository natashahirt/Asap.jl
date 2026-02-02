#=
P-Delta Analysis
================

Iterative analysis that accounts for the second-order effect of gravity loads
acting through lateral displacements.

Algorithm:
1. Solve linear elastic: K·u₁ = F
2. Compute geometric stiffness Kg from internal forces
3. Solve: (K + Kg)·u₂ = F  
4. Update Kg from new forces
5. Repeat until convergence: ||u_{n+1} - u_n|| / ||u_n|| < tolerance

Physical Interpretation:
- Gravity loads "P" acting through lateral displacement "Δ" create 
  additional overturning moment P·Δ
- This reduces effective lateral stiffness
- For stable structures: converges to amplified displacement
- For unstable structures: diverges (buckling)

Code Requirements:
- AISC 360-16: Required when B₂ > 1.1 or second-order effects > 5%
- ASCE 7 Seismic: Stability coefficient θ check

References:
- AISC 360-16 Chapter C and Appendix 7
- McGuire, Gallagher, Ziemian "Matrix Structural Analysis" Ch. 10
=#

using LinearAlgebra: norm
using SparseArrays: SparseMatrixCSC

# =============================================================================
# P-Delta Solver
# =============================================================================

"""
    solve_pdelta!(model; max_iter=10, tol=1e-3, verbose=false) -> PDeltaResult

Perform P-delta analysis with iterative geometric stiffness update.

# Arguments
- `model`: Structural model (FrameModel or Model)
- `max_iter::Int`: Maximum iterations (default: 10)
- `tol::Float64`: Convergence tolerance on displacement (default: 1e-3)
- `verbose::Bool`: Print iteration info (default: false)

# Returns
- `PDeltaResult`: Contains convergence info, iterations, amplification factor

# Algorithm
1. Linear static analysis: K·u₁ = F
2. Compute Kg from internal forces  
3. Solve (K + Kg)·u = F
4. Update internal forces and Kg
5. Repeat until ||Δu||/||u|| < tol

# Example
```julia
model = Model(nodes, columns ∪ beams, gravity_loads ∪ lateral_loads)
result = solve_pdelta!(model; max_iter=10, tol=0.001)

if result.converged
    println("Converged in ", result.iterations, " iterations")
    println("Amplification factor: ", result.amplification)
else
    println("WARNING: Did not converge - structure may be unstable")
end
```

# Notes
- Final displacements and reactions are stored in nodes after analysis
- If diverging, check if structure is near buckling (solve_buckling!)
- Amplification factor B₂ ≈ 1/(1 - P·Δ/(H·V)) where H=height, V=base shear
"""
function solve_pdelta!(model::FrameModel; 
                       max_iter::Int = 10, 
                       tol::Float64 = 1e-3,
                       verbose::Bool = false)
    
    # Initial linear solve
    if !model.processed
        process!(model)
    end
    solve!(model)
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        return PDeltaResult(true, 0, 1.0, 0.0)
    end
    
    # Store first-order displacement for amplification calculation
    u_first = copy(model.u)
    u_prev = copy(model.u)
    
    # Get load vector (constant through iterations)
    F = model.P[idx] - model.Pf[idx]
    
    converged = false
    iter = 0
    
    for i in 1:max_iter
        iter = i
        
        # Assemble geometric stiffness from current internal forces
        Kg = assemble_geometric_stiffness(model)
        
        # Total tangent stiffness
        K_total = model.S + Kg
        
        # Solve with updated stiffness
        U = K_total[idx, idx] \ F
        
        # Update model displacement vector
        model.u = zeros(model.nDOFs)
        model.u[idx] = U
        
        # Post-process to get new internal forces
        post_process!(model)
        
        # Check convergence
        du = model.u - u_prev
        rel_change = norm(du) / max(norm(model.u), 1e-10)
        
        if verbose
            println("Iteration $i: rel_change = $(round(rel_change, sigdigits=4))")
        end
        
        if rel_change < tol
            converged = true
            break
        end
        
        # Check for divergence
        if norm(model.u) > 100 * norm(u_first) && norm(u_first) > 1e-10
            @warn "P-delta iteration diverging - structure may be unstable"
            break
        end
        
        u_prev = copy(model.u)
    end
    
    # Calculate amplification factor
    u_max_first = maximum(abs.(u_first))
    u_max_final = maximum(abs.(model.u))
    amplification = u_max_first > 1e-10 ? u_max_final / u_max_first : 1.0
    
    # Calculate maximum drift ratio (simplified - assumes vertical columns)
    max_drift_ratio = compute_max_drift_ratio(model)
    
    if verbose
        if converged
            println("Converged in $iter iterations")
        else
            println("Did not converge after $max_iter iterations")
        end
        println("Amplification factor: $(round(amplification, digits=4))")
    end
    
    return PDeltaResult(converged, iter, amplification, max_drift_ratio)
end

function solve_pdelta!(model::Model; 
                       max_iter::Int = 10, 
                       tol::Float64 = 1e-3,
                       verbose::Bool = false)
    
    if !model.processed
        process!(model)
    end
    solve!(model)
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        return PDeltaResult(true, 0, 1.0, 0.0)
    end
    
    u_first = copy(model.u)
    u_prev = copy(model.u)
    
    F = model.P[idx] - model.Pf[idx]
    
    converged = false
    iter = 0
    
    for i in 1:max_iter
        iter = i
        
        Kg = assemble_geometric_stiffness(model)
        K_total = model.S + Kg
        
        U = K_total[idx, idx] \ F
        
        model.u = zeros(model.nDOFs)
        model.u[idx] = U
        
        post_process!(model)
        
        du = model.u - u_prev
        rel_change = norm(du) / max(norm(model.u), 1e-10)
        
        if verbose
            println("Iteration $i: rel_change = $(round(rel_change, sigdigits=4))")
        end
        
        if rel_change < tol
            converged = true
            break
        end
        
        if norm(model.u) > 100 * norm(u_first) && norm(u_first) > 1e-10
            @warn "P-delta iteration diverging - structure may be unstable"
            break
        end
        
        u_prev = copy(model.u)
    end
    
    u_max_first = maximum(abs.(u_first))
    u_max_final = maximum(abs.(model.u))
    amplification = u_max_first > 1e-10 ? u_max_final / u_max_first : 1.0
    
    max_drift_ratio = compute_max_drift_ratio(model)
    
    if verbose
        if converged
            println("Converged in $iter iterations")
        else
            println("Did not converge after $max_iter iterations")
        end
        println("Amplification factor: $(round(amplification, digits=4))")
    end
    
    return PDeltaResult(converged, iter, amplification, max_drift_ratio)
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
Compute maximum inter-story drift ratio (Δ/h) - simplified estimate.
"""
function compute_max_drift_ratio(model::AbstractModel)
    max_drift = 0.0
    
    # Get all node positions and displacements
    for node in model.nodes
        # Lateral displacement (X and Y) - strip units
        dx = abs(ustrip(node.displacement[1]))
        dy = abs(ustrip(node.displacement[2]))
        
        # Height (Z coordinate)
        z = ustrip(u"m", node.position[3])
        
        if z > 0.1  # Avoid division by very small heights
            drift_x = dx / z
            drift_y = dy / z
            max_drift = max(max_drift, drift_x, drift_y)
        end
    end
    
    return max_drift
end

"""
    amplification_factor(model) -> Float64

Convenience function to get just the P-delta amplification factor.
"""
function amplification_factor(model::Union{FrameModel, Model})
    result = solve_pdelta!(model)
    return result.amplification
end

"""
    B2_factor(model) -> Float64

Calculate the AISC B₂ amplification factor from P-delta analysis.
Equivalent to amplification_factor but named for code familiarity.
"""
B2_factor(model) = amplification_factor(model)
