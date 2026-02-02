#=
Nonlinear Static Analysis (Newton-Raphson)
==========================================

Full nonlinear static analysis with load increments and Newton-Raphson iteration.

Algorithm (Load-Controlled):
For each load step λ:
    1. Apply load increment: F = λ·F_total
    2. Newton-Raphson iteration until equilibrium:
       a. Compute residual: R = F_ext - F_int(u)
       b. Compute tangent stiffness: Kt = K + Kg
       c. Solve for increment: Kt·Δu = R
       d. Update: u = u + Δu
       e. Check convergence: ||R|| < tol
    3. Store converged state
    4. Continue to next load step

Applications:
- Pushover analysis for seismic assessment
- Capacity analysis (load-deflection curves)
- Structures with large displacements
- Post-buckling behavior (with arc-length methods)

References:
- Cook et al. "Concepts and Applications of FEA" Ch. 16
- Crisfield "Non-linear Finite Element Analysis of Solids and Structures"
- ASCE 41 - Seismic Evaluation and Retrofit (nonlinear static procedure)
=#

using LinearAlgebra: norm
using SparseArrays: SparseMatrixCSC

# =============================================================================
# Nonlinear Static Solver
# =============================================================================

"""
    solve_nonlinear!(model; n_steps=10, max_iter=20, tol=1e-4, verbose=false) -> NonlinearResult

Perform nonlinear static analysis with incremental loading and Newton-Raphson iteration.

# Arguments
- `model`: Structural model (FrameModel or Model)
- `n_steps::Int`: Number of load increments (default: 10)
- `max_iter::Int`: Max Newton-Raphson iterations per step (default: 20)
- `tol::Float64`: Equilibrium tolerance (default: 1e-4)
- `verbose::Bool`: Print iteration info (default: false)

# Returns
- `NonlinearResult`: Load factors, displacements at each step, iteration counts

# Example
```julia
# Pushover analysis
model = Model(nodes, elements, [lateral_load])
result = solve_nonlinear!(model; n_steps=20, tol=1e-5)

# Plot capacity curve
using Plots
plot([maximum(abs.(u)) for u in result.displacements], result.load_factors,
     xlabel="Max Displacement", ylabel="Load Factor")
```

# Notes
- Uses load-controlled increments (not arc-length)
- May fail to converge near limit points (use finer steps)
- Final state is at full load (λ = 1.0)
- Intermediate results stored for plotting load-deflection curves
"""
function solve_nonlinear!(model::FrameModel; 
                          n_steps::Int = 10,
                          max_iter::Int = 20,
                          tol::Float64 = 1e-4,
                          verbose::Bool = false)
    
    # Process model
    if !model.processed
        process!(model)
    end
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        return NonlinearResult(true, [1.0], [zeros(model.nDOFs)], [0], [0.0])
    end
    
    # Total external load (what we're incrementing toward)
    F_total = model.P[idx] - model.Pf[idx]
    
    # Initialize
    model.u = zeros(model.nDOFs)
    
    # Storage for results
    load_factors = Float64[]
    displacements = Vector{Float64}[]
    iterations_per_step = Int[]
    equilibrium_errors = Float64[]
    
    # Initial state
    push!(load_factors, 0.0)
    push!(displacements, copy(model.u))
    push!(iterations_per_step, 0)
    push!(equilibrium_errors, 0.0)
    
    # Load increments
    dλ = 1.0 / n_steps
    λ = 0.0
    
    overall_converged = true
    
    for step in 1:n_steps
        λ += dλ
        F_ext = λ * F_total
        
        if verbose
            println("\n--- Load Step $step: λ = $(round(λ, digits=4)) ---")
        end
        
        # Newton-Raphson iteration
        converged = false
        iter_count = 0
        residual_norm = 0.0
        
        for iter in 1:max_iter
            iter_count = iter
            
            # Compute internal forces from current displacement
            # (This updates element.forces via post_process!)
            post_process!(model)
            
            # Assemble tangent stiffness: Kt = K + Kg
            Kg = assemble_geometric_stiffness(model)
            Kt = model.S + Kg
            
            # Compute internal force vector
            F_int = compute_internal_forces(model, idx)
            
            # Residual (out-of-balance force)
            R = F_ext - F_int
            residual_norm = norm(R)
            
            if verbose
                println("  Iter $iter: ||R|| = $(round(residual_norm, sigdigits=4))")
            end
            
            # Check convergence
            if residual_norm < tol * max(norm(F_ext), 1.0)
                converged = true
                break
            end
            
            # Solve for displacement increment
            Kt_ff = Kt[idx, idx]
            Δu = Kt_ff \ R
            
            # Update displacement
            model.u[idx] += Δu
            
            # Check for divergence
            if norm(Δu) > 1e6 || any(isnan.(Δu))
                @warn "Newton-Raphson diverging at step $step, iteration $iter"
                break
            end
        end
        
        if !converged
            @warn "Did not converge at load step $step (λ = $λ)"
            overall_converged = false
        end
        
        # Store results for this step
        push!(load_factors, λ)
        push!(displacements, copy(model.u))
        push!(iterations_per_step, iter_count)
        push!(equilibrium_errors, residual_norm)
        
        # Update node displacements for visualization
        post_process!(model)
    end
    
    if verbose
        println("\n=== Nonlinear Analysis Complete ===")
        println("Converged: $overall_converged")
        println("Total iterations: $(sum(iterations_per_step))")
    end
    
    return NonlinearResult(overall_converged, load_factors, displacements, 
                           iterations_per_step, equilibrium_errors)
end

function solve_nonlinear!(model::Model; 
                          n_steps::Int = 10,
                          max_iter::Int = 20,
                          tol::Float64 = 1e-4,
                          verbose::Bool = false)
    
    if !model.processed
        process!(model)
    end
    
    idx = model.freeDOFs
    n_free = length(idx)
    
    if n_free == 0
        return NonlinearResult(true, [1.0], [zeros(model.nDOFs)], [0], [0.0])
    end
    
    F_total = model.P[idx] - model.Pf[idx]
    
    model.u = zeros(model.nDOFs)
    
    load_factors = Float64[]
    displacements = Vector{Float64}[]
    iterations_per_step = Int[]
    equilibrium_errors = Float64[]
    
    push!(load_factors, 0.0)
    push!(displacements, copy(model.u))
    push!(iterations_per_step, 0)
    push!(equilibrium_errors, 0.0)
    
    dλ = 1.0 / n_steps
    λ = 0.0
    
    overall_converged = true
    
    for step in 1:n_steps
        λ += dλ
        F_ext = λ * F_total
        
        if verbose
            println("\n--- Load Step $step: λ = $(round(λ, digits=4)) ---")
        end
        
        converged = false
        iter_count = 0
        residual_norm = 0.0
        
        for iter in 1:max_iter
            iter_count = iter
            
            post_process!(model)
            
            Kg = assemble_geometric_stiffness(model)
            Kt = model.S + Kg
            
            F_int = compute_internal_forces(model, idx)
            
            R = F_ext - F_int
            residual_norm = norm(R)
            
            if verbose
                println("  Iter $iter: ||R|| = $(round(residual_norm, sigdigits=4))")
            end
            
            if residual_norm < tol * max(norm(F_ext), 1.0)
                converged = true
                break
            end
            
            Kt_ff = Kt[idx, idx]
            Δu = Kt_ff \ R
            
            model.u[idx] += Δu
            
            if norm(Δu) > 1e6 || any(isnan.(Δu))
                @warn "Newton-Raphson diverging at step $step, iteration $iter"
                break
            end
        end
        
        if !converged
            @warn "Did not converge at load step $step (λ = $λ)"
            overall_converged = false
        end
        
        push!(load_factors, λ)
        push!(displacements, copy(model.u))
        push!(iterations_per_step, iter_count)
        push!(equilibrium_errors, residual_norm)
        
        post_process!(model)
    end
    
    if verbose
        println("\n=== Nonlinear Analysis Complete ===")
        println("Converged: $overall_converged")
        println("Total iterations: $(sum(iterations_per_step))")
    end
    
    return NonlinearResult(overall_converged, load_factors, displacements, 
                           iterations_per_step, equilibrium_errors)
end

# =============================================================================
# Internal Force Computation
# =============================================================================

"""
Compute the internal force vector from current element forces.
"""
function compute_internal_forces(model::FrameModel, free_dofs::Vector{Int})
    F_int = zeros(model.nDOFs)
    
    for element in model.elements
        if !(element isa FrameElement)
            continue
        end
        
        # Element forces in global coordinates
        # F_global = R' * F_local
        R = element.R
        f_local = element.forces
        f_global = R' * f_local
        
        # Assemble into global vector
        gid = element.globalID
        for i in 1:12
            F_int[gid[i]] += f_global[i]
        end
    end
    
    return F_int[free_dofs]
end

function compute_internal_forces(model::Model, free_dofs::Vector{Int})
    F_int = zeros(model.nDOFs)
    
    # Frame elements
    for element in model.frame_elements
        if !(element isa FrameElement)
            continue
        end
        
        R = element.R
        f_local = element.forces
        f_global = R' * f_local
        
        gid = element.globalID
        for i in 1:12
            F_int[gid[i]] += f_global[i]
        end
    end
    
    # Shell elements - for now, use K*u approximation
    # (Full shell internal force requires stress integration)
    for shell in model.shell_elements
        # Shell contribution: F_int = K_shell * u_shell
        gid = shell.globalID
        u_elem = model.u[gid]
        f_elem = shell.K * u_elem
        for i in eachindex(gid)
            F_int[gid[i]] += f_elem[i]
        end
    end
    
    return F_int[free_dofs]
end

# =============================================================================
# Convenience Functions
# =============================================================================

"""
    pushover!(model; kwargs...) -> NonlinearResult

Alias for solve_nonlinear! - clearer name for seismic assessment context.
"""
pushover!(model; kwargs...) = solve_nonlinear!(model; kwargs...)

"""
    capacity_curve(result::NonlinearResult) -> NamedTuple{(:displacement, :load_factor)}

Extract the load-displacement capacity curve from nonlinear analysis.

Returns (displacement=..., load_factor=...) suitable for plotting.

# Example
```julia
result = solve!(model, :nonlinear)
curve = capacity_curve(result)
plot(curve.displacement, curve.load_factor)
```
"""
function capacity_curve(result::NonlinearResult)
    max_disps = [maximum(abs.(u)) for u in result.displacements]
    return (displacement=max_disps, load_factor=result.load_factors)
end
