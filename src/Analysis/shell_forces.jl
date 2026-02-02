#=
Shell Internal Forces
=====================
Internal force recovery for shell elements (ShellTri3, CompositeShellTri3).

Attribution: Shell element formulation based on FinEtoolsFlexStructures.jl
by Petr Krysl (UC San Diego).

Shell resultants are given per unit length at the element centroid:
- Membrane forces: Nxx, Nyy, Nxy [N/m]
- Bending moments: Mxx, Myy, Mxy [N·m/m = N]  
- Transverse shear: Qxz, Qyz [N/m]

These are stress resultants integrated through the thickness.
=#

# Formatting helper
_fmt(x::Float64, decimals::Int=2) = string(round(x, digits=decimals))

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

# ============================================================================
# ShellInternalForces
# ============================================================================

"""
    ShellInternalForces

Internal force resultants for a shell element evaluated at its centroid.

# Fields
- `element::ShellElement` - Reference to the shell element
- `Nxx::Float64` - Membrane force in local X direction [N/m]
- `Nyy::Float64` - Membrane force in local Y direction [N/m]
- `Nxy::Float64` - Membrane shear force [N/m]
- `Mxx::Float64` - Bending moment about local Y axis [N·m/m]
- `Myy::Float64` - Bending moment about local X axis [N·m/m]
- `Mxy::Float64` - Twisting moment [N·m/m]
- `Qxz::Float64` - Transverse shear in X-Z plane [N/m]
- `Qyz::Float64` - Transverse shear in Y-Z plane [N/m]
- `ply_stresses::Union{Nothing, Vector{Vector{Float64}}}` - Ply stresses [σ1, σ2, τ12] 
   in ply material coordinates for each ply (composites only, `nothing` for isotropic)

# Sign Convention
- Positive Mxx causes tension on +Z face
- Positive Qxz acts in +Z direction on +X face
"""
struct ShellInternalForces
    element::ShellElement
    Nxx::Float64
    Nyy::Float64
    Nxy::Float64
    Mxx::Float64
    Myy::Float64
    Mxy::Float64
    Qxz::Float64
    Qyz::Float64
    ply_stresses::Union{Nothing, Vector{Vector{Float64}}}
end

# Type stability
Base.eltype(::Type{Vector{ShellInternalForces}}) = ShellInternalForces

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
    u_local = _extract_local_dofs(elem, u_global)
    
    # Compute gradient matrix and area
    gradN_e = zeros(3, 2)
    gradN_e, Ae = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    # --- Membrane forces ---
    Bm = zeros(3, 18)
    _Bmmat!(Bm, gradN_e)
    ε_membrane = Bm * u_local
    
    Dps = _plane_stress_D(elem.E, elem.ν)
    σ_membrane = Dps * ε_membrane
    
    t = elem.thickness
    Nxx = σ_membrane[1] * t
    Nyy = σ_membrane[2] * t
    Nxy = σ_membrane[3] * t
    
    # --- Bending moments ---
    Bb = zeros(3, 18)
    _Bbmat!(Bb, gradN_e)
    κ = Bb * u_local
    
    D_bend = (t^3 / 12.0) * Dps
    M = D_bend * κ
    Mxx, Myy, Mxy = M[1], M[2], M[3]
    
    # --- Transverse shear ---
    Bs = zeros(2, 18)
    _Bsmat!(Bs, elem.ecoords_e, Ae)
    γ = Bs * u_local
    
    G = elem.E / (2 * (1 + elem.ν))
    κ_shear = 5.0 / 6.0
    
    Qxz = κ_shear * G * t * γ[1]
    Qyz = κ_shear * G * t * γ[2]
    
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
    [ShellInternalForces(elem, model.u) for elem in model.shells]
end

function ShellInternalForces(model::Model)
    [ShellInternalForces(elem, model.u) for elem in model.shells]
end

# ============================================================================
# Principal Values and Design Quantities
# ============================================================================

"""
    principal_moments(sif::ShellInternalForces)

Compute principal bending moments and their directions.

# Returns
- `(M1=..., M2=..., θ=...)`: NamedTuple with principal moments and angle (radians)
"""
function principal_moments(sif::ShellInternalForces)
    Mx, My, Mxy = sif.Mxx, sif.Myy, sif.Mxy
    
    M_avg = (Mx + My) / 2
    R = sqrt(((Mx - My) / 2)^2 + Mxy^2)
    
    M1 = M_avg + R
    M2 = M_avg - R
    θ = 0.5 * atan(2 * Mxy, Mx - My)
    
    return (M1=M1, M2=M2, θ=θ)
end

"""
    principal_forces(sif::ShellInternalForces)

Compute principal membrane forces and their directions.

# Returns
- `(N1=..., N2=..., θ=...)`: NamedTuple with principal forces and angle (radians)
"""
function principal_forces(sif::ShellInternalForces)
    Nx, Ny, Nxy = sif.Nxx, sif.Nyy, sif.Nxy
    
    N_avg = (Nx + Ny) / 2
    R = sqrt(((Nx - Ny) / 2)^2 + Nxy^2)
    
    N1 = N_avg + R
    N2 = N_avg - R
    θ = 0.5 * atan(2 * Nxy, Nx - Ny)
    
    return (N1=N1, N2=N2, θ=θ)
end

"""
    von_mises_stress(sif::ShellInternalForces, z::Float64, t::Float64)

Compute von Mises stress at a given z-coordinate through the thickness.

# Arguments
- `z`: Distance from mid-plane (positive = top face)
- `t`: Total thickness

# Returns
von Mises stress in Pa
"""
function von_mises_stress(sif::ShellInternalForces, z::Float64, t::Float64)
    σxx = sif.Nxx / t + 12 * sif.Mxx * z / t^3
    σyy = sif.Nyy / t + 12 * sif.Myy * z / t^3
    τxy = sif.Nxy / t + 12 * sif.Mxy * z / t^3
    
    sqrt(σxx^2 - σxx*σyy + σyy^2 + 3*τxy^2)
end

"""
    max_surface_stresses(sif::ShellInternalForces, t::Float64)

Compute maximum stresses on top (+t/2) and bottom (-t/2) surfaces.

# Returns
Named tuple: `(top=(σxx, σyy, τxy), bottom=(σxx, σyy, τxy))`
"""
function max_surface_stresses(sif::ShellInternalForces, t::Float64)
    σxx_top = sif.Nxx / t + 6 * sif.Mxx / t^2
    σyy_top = sif.Nyy / t + 6 * sif.Myy / t^2
    τxy_top = sif.Nxy / t + 6 * sif.Mxy / t^2
    
    σxx_bot = sif.Nxx / t - 6 * sif.Mxx / t^2
    σyy_bot = sif.Nyy / t - 6 * sif.Myy / t^2
    τxy_bot = sif.Nxy / t - 6 * sif.Mxy / t^2
    
    (top = (σxx=σxx_top, σyy=σyy_top, τxy=τxy_top),
     bottom = (σxx=σxx_bot, σyy=σyy_bot, τxy=τxy_bot))
end

# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, sif::ShellInternalForces)
    print(io, "ShellInternalForces(")
    print(io, "N=[$(_fmt(sif.Nxx, 1)), $(_fmt(sif.Nyy, 1)), $(_fmt(sif.Nxy, 1))], ")
    print(io, "M=[$(_fmt(sif.Mxx, 2)), $(_fmt(sif.Myy, 2)), $(_fmt(sif.Mxy, 2))], ")
    print(io, "Q=[$(_fmt(sif.Qxz, 1)), $(_fmt(sif.Qyz, 1))])")
end

function Base.show(io::IO, ::MIME"text/plain", sif::ShellInternalForces)
    id = hasproperty(sif.element, :id) ? sif.element.id : "unknown"
    println(io, "ShellInternalForces for element ':$(id)'")
    println(io, "├─ Membrane forces [N/m]:")
    println(io, "│    Nxx = $(_fmt(sif.Nxx, 2))")
    println(io, "│    Nyy = $(_fmt(sif.Nyy, 2))")
    println(io, "│    Nxy = $(_fmt(sif.Nxy, 2))")
    println(io, "├─ Bending moments [N·m/m]:")
    println(io, "│    Mxx = $(_fmt(sif.Mxx, 4))")
    println(io, "│    Myy = $(_fmt(sif.Myy, 4))")
    println(io, "│    Mxy = $(_fmt(sif.Mxy, 4))")
    println(io, "└─ Transverse shear [N/m]:")
    println(io, "     Qxz = $(_fmt(sif.Qxz, 2))")
    print(io, "     Qyz = $(_fmt(sif.Qyz, 2))")
end

# ============================================================================
# ShellDisplacements
# ============================================================================

"""
    ShellDisplacements

Displacement results for a shell element at its centroid.

# Fields
- `element::ShellElement` - Reference to the shell element
- `u::Float64` - Translation in local X [m]
- `v::Float64` - Translation in local Y [m]
- `w::Float64` - Translation in local Z (out-of-plane) [m]
- `θx::Float64` - Rotation about local X [rad]
- `θy::Float64` - Rotation about local Y [rad]
- `θz::Float64` - Rotation about local Z (drilling) [rad]
"""
struct ShellDisplacements
    element::ShellElement
    u::Float64
    v::Float64
    w::Float64
    θx::Float64
    θy::Float64
    θz::Float64
end

# Type stability
Base.eltype(::Type{Vector{ShellDisplacements}}) = ShellDisplacements

"""
    ShellDisplacements(elem::ShellElement, u_global::Vector{Float64})

Compute shell displacement at element centroid.
Works for all shell element types (ShellTri3, CompositeShellTri3).
"""
function ShellDisplacements(elem::ShellElement, u_global::Vector{Float64})
    u_local = _extract_local_dofs(elem, u_global)
    
    # Centroid at ξ = η = 1/3: N1 = N2 = N3 = 1/3
    N = 1.0 / 3.0
    
    u = N * (u_local[1] + u_local[7] + u_local[13])
    v = N * (u_local[2] + u_local[8] + u_local[14])
    w = N * (u_local[3] + u_local[9] + u_local[15])
    θx = N * (u_local[4] + u_local[10] + u_local[16])
    θy = N * (u_local[5] + u_local[11] + u_local[17])
    θz = N * (u_local[6] + u_local[12] + u_local[18])
    
    ShellDisplacements(elem, u, v, w, θx, θy, θz)
end

"""Compute displacements for all shell elements in a model."""
function ShellDisplacements(model::ShellModel)
    [ShellDisplacements(elem, model.u) for elem in model.shells]
end

function ShellDisplacements(model::Model)
    [ShellDisplacements(elem, model.u) for elem in model.shells]
end

"""
    max_deflection(disps::Vector{ShellDisplacements})

Find the maximum out-of-plane deflection across all shells.

# Returns
Named tuple: `(max_w, element)` with max deflection value and corresponding element.
"""
function max_deflection(disps::Vector{ShellDisplacements})
    isempty(disps) && return (max_w=0.0, element=nothing)
    max_idx = argmax(abs.(getfield.(disps, :w)))
    (max_w=disps[max_idx].w, element=disps[max_idx].element)
end

function Base.show(io::IO, sd::ShellDisplacements)
    id = hasproperty(sd.element, :id) ? sd.element.id : "unknown"
    print(io, "ShellDisplacements(:$id, w=$(_fmt(sd.w*1000, 3))mm)")
end

function Base.show(io::IO, ::MIME"text/plain", sd::ShellDisplacements)
    id = hasproperty(sd.element, :id) ? sd.element.id : "unknown"
    println(io, "ShellDisplacements for element ':$(id)'")
    println(io, "├─ Translations [m]:")
    println(io, "│    u = $(_fmt(sd.u, 6))")
    println(io, "│    v = $(_fmt(sd.v, 6))")
    println(io, "│    w = $(_fmt(sd.w, 6))")
    println(io, "└─ Rotations [rad]:")
    println(io, "     θx = $(_fmt(sd.θx, 6))")
    println(io, "     θy = $(_fmt(sd.θy, 6))")
    print(io, "     θz = $(_fmt(sd.θz, 6))")
end

# ============================================================================
# Unified InternalForces Interface
# ============================================================================

"""
    InternalForces(element, model; kwargs...)

Unified interface for internal force recovery. Dispatches to:
- `ElementInternalForces` for frame elements (beams, columns, trusses)
- `ShellInternalForces` for shell elements

# Example
```julia
model = Model(nodes, frames, shells, [load])
process!(model)
solve!(model)

beam_forces = InternalForces(frames[1], model)   # → ElementInternalForces
shell_forces = InternalForces(shells[1], model)  # → ShellInternalForces
```
"""
function InternalForces(elem::FrameElement, model::AbstractModel; resolution=20)
    ElementInternalForces(elem, model; resolution=resolution)
end

function InternalForces(elem::ShellElement, model::AbstractModel)
    ShellInternalForces(elem, model.u)
end

"""
    InternalForces(model; increment=0.1)

Get internal forces for all elements in a model.

# Returns
Named tuple: `(frames=Vector{ElementInternalForces}, shells=Vector{ShellInternalForces})`
"""
function InternalForces(model::Model; increment=0.1)
    frame_forces = ElementInternalForces(model, increment)
    shell_forces = ShellInternalForces(model)
    (frames=frame_forces, shells=shell_forces)
end

function InternalForces(model::FrameModel; increment=0.1)
    (frames=ElementInternalForces(model, increment), shells=ShellInternalForces[])
end

function InternalForces(model::ShellModel)
    (frames=ElementInternalForces[], shells=ShellInternalForces(model))
end

# ============================================================================
# Unified Displacements Interface
# ============================================================================

"""
    Displacements(element, model)

Get displacement for a single element.

# Returns
- `ShellDisplacements` for shell elements
- For frame elements, use `displacements(element, model)` instead

# Example
```julia
disp = Displacements(shell1, model)
println("w = \$(disp.w) m")
```
"""
function Displacements(elem::ShellElement, model::AbstractModel)
    ShellDisplacements(elem, model.u)
end

"""
    Displacements(model)

Get displacements for all elements in a model.

# Returns
Named tuple: `(shells=Vector{ShellDisplacements},)`

# Example
```julia
disps = Displacements(model)
disps.shells  # All shell centroid displacements
```
"""
function Displacements(model::Model)
    (shells=ShellDisplacements(model),)
end

function Displacements(model::ShellModel)
    (shells=ShellDisplacements(model),)
end
