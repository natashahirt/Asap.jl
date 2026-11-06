# ============================================================================
# Helper: Extract element DOFs from global displacement vector
# ============================================================================

"""
    _extract_local_dofs(elem::ShellElement, u_global::Vector{Float64})

Extract and transform element DOFs from global displacement vector to local coordinates.
Returns the 18-element local displacement vector.
"""
function _extract_local_dofs(elem::ShellElement, u_global::Vector{Float64})
    u_elem = zeros(18)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:6
            u_elem[(i-1)*6 + j] = u_global[node.globalID[j]]
        end
    end
    return elem.R * u_elem
end

"""In-place variant: fills `u_elem` and `u_local` from workspace."""
@inline function _extract_local_dofs!(u_elem::Vector{Float64}, u_local::Vector{Float64},
                                       elem::ShellElement, u_global::Vector{Float64})
    fill!(u_elem, 0.0)
    @inbounds for (i, node) in enumerate(elem.nodes)
        off = (i - 1) * 6
        for j in 1:6
            u_elem[off + j] = u_global[node.globalID[j]]
        end
    end
    mul!(u_local, elem.R, u_elem)
    return u_local
end

"""
    ShellInternalForces(elem::ShellTri3, u_global::Vector{Float64})

Compute all shell internal force resultants for a ShellTri3 element.

# Arguments
- `elem::ShellTri3` - Processed shell element (must have R, globalID populated)
- `u_global::Vector{Float64}` - Global displacement vector from solved model

# Returns
- `ShellInternalForces` - Internal force resultants at element centroid

# Example
```julia
model = Model(nodes, [], shells, [load])
process!(model)
solve!(model)

sif = ShellInternalForces(shells[1], model.u)
println("Mx = \$(sif.Mxx) N·m/m")
```
"""
function ShellInternalForces(elem::ShellTri3, u_global::Vector{Float64})
    ws = _thread_buffer(_SHELL_FORCES_WS)
    return ShellInternalForces(elem, u_global, ws)
end

"""
    ShellInternalForces(elem::ShellTri3, u_global, ws::ShellForcesWorkspace)

Zero-allocation constructor using pre-allocated workspace buffers.
"""
function ShellInternalForces(elem::ShellTri3, u_global::Vector{Float64}, ws::ShellForcesWorkspace)
    _extract_local_dofs!(ws.u_elem, ws.u_local, elem, u_global)
    
    # Compute gradient matrix and area
    _, Ae = _compute_gradN_Ae!(ws.gradN_e, elem.ecoords_e)
    
    # --- Membrane forces ---
    _Bmmat!(ws.Bm, ws.gradN_e)        # _Bmmat! calls fill! internally
    mul!(ws.ε, ws.Bm, ws.u_local)
    
    _plane_stress_D!(ws.Dps, elem.E, elem.ν)
    mul!(ws.σ, ws.Dps, ws.ε)
    
    t = elem.thickness
    Nxx = ws.σ[1] * t
    Nyy = ws.σ[2] * t
    Nxy = ws.σ[3] * t
    
    # --- Bending moments ---
    _Bbmat!(ws.Bb, ws.gradN_e)         # _Bbmat! calls fill! internally
    mul!(ws.κ, ws.Bb, ws.u_local)
    
    fac = t^3 / 12.0
    @inbounds for j in 1:3; ws.M[j] = 0.0; end
    @inbounds for j in 1:3, i in 1:3
        ws.M[i] += fac * ws.Dps[i, j] * ws.κ[j]
    end
    Mxx, Myy, Mxy = ws.M[1], ws.M[2], ws.M[3]
    
    # --- Transverse shear ---
    _Bsmat!(ws.Bs, elem.ecoords_e, Ae) # _Bsmat! fills internally
    mul!(ws.γ, ws.Bs, ws.u_local)
    
    G = elem.E / (2 * (1 + elem.ν))
    κ_shear = 5.0 / 6.0
    
    Qxz = κ_shear * G * t * ws.γ[1]
    Qyz = κ_shear * G * t * ws.γ[2]
    
    return ShellInternalForces(elem, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qxz, Qyz, nothing)
end

"""
    ShellInternalForces(elem::CompositeShellTri3, u_global::Vector{Float64})

Compute internal force resultants for a composite shell element.
Uses the laminate ABD matrix for proper coupling between membrane and bending.
Also computes ply stresses at the middle of each ply.
"""
function ShellInternalForces(elem::CompositeShellTri3, u_global::Vector{Float64})
    u_local = _extract_local_dofs(elem, u_global)
    
    # Compute gradient matrix and area
    gradN_e = zeros(3, 2)
    gradN_e, Ae = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    # --- Membrane strains ---
    Bm = zeros(3, 18)
    _Bmmat!(Bm, gradN_e)
    ε0 = Bm * u_local
    
    # --- Curvatures ---
    Bb = zeros(3, 18)
    _Bbmat!(Bb, gradN_e)
    κ = Bb * u_local
    
    # Get laminate ABD matrix
    A, B, D = elem.laminate.A, elem.laminate.B, elem.laminate.D
    
    # Resultant forces and moments
    N = A * ε0 + B * κ
    M = B * ε0 + D * κ
    
    Nxx, Nyy, Nxy = N[1], N[2], N[3]
    Mxx, Myy, Mxy = M[1], M[2], M[3]
    
    # --- Transverse shear (effective properties) ---
    Bs = zeros(2, 18)
    _Bsmat!(Bs, elem.ecoords_e, Ae)
    γ = Bs * u_local
    
    h = elem.laminate.h
    G_eff = elem.laminate.E_eq / (2 * (1 + elem.laminate.ν_eq))
    κ_shear = 5.0 / 6.0
    
    Qxz = κ_shear * G_eff * h * γ[1]
    Qyz = κ_shear * G_eff * h * γ[2]
    
    # --- Compute ply stresses at middle of each ply ---
    n_plies = length(elem.laminate.plies)
    ply_stress_results = Vector{Vector{Float64}}(undef, n_plies)
    
    for i in 1:n_plies
        # Use the existing ply_stresses function from composite.jl
        ply_stress_results[i] = ply_stresses(elem, u_global, i; z_position=:middle)
    end
    
    return ShellInternalForces(elem, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qxz, Qyz, ply_stress_results)
end

"""
    ShellInternalForces(model::AbstractModel)

Compute internal forces for all shell elements in a model.

# Returns
Vector of `ShellInternalForces`, one per shell element.
"""
function ShellInternalForces(model::ShellModel)
    [ShellInternalForces(elem, model.u) for elem in model.elements]
end

function ShellInternalForces(model::Model)
    [ShellInternalForces(elem, model.u) for elem in model.shell_elements]
end

