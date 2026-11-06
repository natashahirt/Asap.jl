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
# ShellForcesWorkspace — pre-allocated buffers for zero-alloc ShellInternalForces
# ============================================================================

"""
    ShellForcesWorkspace

Pre-allocated buffers for computing `ShellInternalForces` without heap allocations.
Create one per thread and pass to `ShellInternalForces(elem, u, ws)`.
"""
struct ShellForcesWorkspace
    u_elem::Vector{Float64}   # 18
    u_local::Vector{Float64}  # 18
    gradN_e::Matrix{Float64}  # 3×2
    Bm::Matrix{Float64}       # 3×18
    Bb::Matrix{Float64}       # 3×18
    Bs::Matrix{Float64}       # 2×18
    ε::Vector{Float64}        # 3
    κ::Vector{Float64}        # 3
    γ::Vector{Float64}        # 2
    Dps::Matrix{Float64}      # 3×3
    σ::Vector{Float64}        # 3
    M::Vector{Float64}        # 3
end

function ShellForcesWorkspace()
    ShellForcesWorkspace(
        zeros(18), zeros(18),
        zeros(3, 2), zeros(3, 18), zeros(3, 18), zeros(2, 18),
        zeros(3), zeros(3), zeros(2),
        zeros(3, 3), zeros(3), zeros(3),
    )
end

# Per-thread pool, filled in `Asap.__init__` (see threading.jl).
const _SHELL_FORCES_WS = ShellForcesWorkspace[]

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

