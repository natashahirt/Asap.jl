#=
Result Types for Nonlinear Analysis
===================================

Container types for analysis results, following the pattern of ModalResult.
=#

using Printf

"""
    BucklingResult

Result container for linear buckling analysis.

# Fields
- `load_factors::Vector{Float64}`: Critical load multipliers (λ)
- `mode_shapes::Matrix{Float64}`: Buckling mode shapes (columns are modes)
- `n_modes::Int`: Number of modes computed

# Interpretation
- `load_factors[1]` is the lowest buckling load multiplier
- If λ₁ > 1.0, structure won't buckle under current loads
- If λ₁ < 1.0, structure will buckle before reaching current load level

# Example
```julia
result = solve_buckling!(model; n=6)
println("Critical load factor: ", result.load_factors[1])
println("Structure is ", result.load_factors[1] > 1.0 ? "stable" : "UNSTABLE")
```
"""
struct BucklingResult
    load_factors::Vector{Float64}    # Critical load multipliers
    mode_shapes::Matrix{Float64}     # Buckling modes (columns)
    n_modes::Int                     # Number of modes computed
end

"""
    PDeltaResult

Result container for P-delta analysis.

# Fields
- `converged::Bool`: Did the iteration converge?
- `iterations::Int`: Number of iterations taken
- `amplification::Float64`: Second-order amplification factor (δ_final / δ_first)
- `max_drift_ratio::Float64`: Maximum story drift ratio

# Example
```julia
result = solve_pdelta!(model; max_iter=10, tol=1e-3)
if result.converged
    println("Amplification factor: ", result.amplification)
end
```
"""
struct PDeltaResult
    converged::Bool
    iterations::Int
    amplification::Float64           # δ₂ / δ₁ overall amplification
    max_drift_ratio::Float64         # Maximum Δ/h
end

"""
    NonlinearResult

Result container for nonlinear static analysis (Newton-Raphson).

# Fields
- `converged::Bool`: Did analysis converge at all load steps?
- `load_factors::Vector{Float64}`: Load multipliers at each step
- `displacements::Vector{Vector{Float64}}`: Displacement vectors at each step
- `iterations_per_step::Vector{Int}`: Iterations needed at each step
- `equilibrium_error::Vector{Float64}`: Final residual norm at each step

# Example
```julia
result = solve_nonlinear!(model; n_steps=10)
# Access final displacement
u_final = result.displacements[end]
# Plot load-displacement curve
plot(result.displacements, result.load_factors)
```
"""
struct NonlinearResult
    converged::Bool
    load_factors::Vector{Float64}           # λ at each step
    displacements::Vector{Vector{Float64}}  # u at each step
    iterations_per_step::Vector{Int}
    equilibrium_error::Vector{Float64}      # ||R|| at convergence
end

# =============================================================================
# Pretty Printing
# =============================================================================

"""Print a formatted summary of buckling analysis results."""
function print_buckling_summary(result::BucklingResult)
    println("Buckling Analysis Summary")
    println("=" ^ 50)
    println("Number of modes: $(result.n_modes)")
    println()
    println("Mode | Load Factor (λ) | Status")
    println("-" ^ 50)
    for i in 1:result.n_modes
        λ = result.load_factors[i]
        status = λ > 1.0 ? "Stable" : "CRITICAL"
        λ_str = lpad(round(λ, digits=4), 14)
        println(" $(lpad(i, 3)) | $λ_str    | $status")
    end
    println()
    if result.load_factors[1] < 1.0
        println("⚠️  WARNING: Structure will buckle before reaching design load!")
    else
        println("✓  Structure is stable under current loading.")
    end
end

"""Print a formatted summary of P-delta analysis results."""
function print_pdelta_summary(result::PDeltaResult)
    println("P-Delta Analysis Summary")
    println("=" ^ 50)
    println("Converged: $(result.converged ? "Yes" : "No")")
    println("Iterations: $(result.iterations)")
    println("Amplification factor: $(round(result.amplification, digits=4))")
    println("Max drift ratio: $(round(result.max_drift_ratio * 100, digits=2))%")
    println()
    if result.amplification > 1.5
        println("⚠️  WARNING: Significant P-delta effects (B₂ > 1.5)")
    elseif result.amplification > 1.1
        println("⚠️  P-delta effects are notable (B₂ > 1.1)")
    else
        println("✓  P-delta effects are minor.")
    end
end

"""Print a formatted summary of nonlinear static analysis results."""
function print_nonlinear_summary(result::NonlinearResult)
    println("Nonlinear Static Analysis Summary")
    println("=" ^ 50)
    println("Converged: $(result.converged ? "Yes" : "No")")
    println("Load steps: $(length(result.load_factors))")
    println("Total iterations: $(sum(result.iterations_per_step))")
    println()
    println("Step | Load Factor | Iterations | Residual")
    println("-" ^ 50)
    for i in eachindex(result.load_factors)
        λ = round(result.load_factors[i], digits=3)
        iter = result.iterations_per_step[i]
        res = @sprintf("%.2e", result.equilibrium_error[i])
        println(" $(lpad(i, 3)) |    $(lpad(λ, 7))  |     $(lpad(iter, 3))    | $res")
    end
end

