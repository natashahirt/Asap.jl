#=================================
Composite Layup Module for Asap.jl
==================================

Core data structures and CLT mathematics for composite laminate analysis.
Provides Ply and Laminate types, plus ABD matrix computation.

Application-specific material presets (carbon, glass, CLT, etc.) should
be defined in your project code, not here.

References:
[1] R.M. Jones, "Mechanics of Composite Materials", 2nd Ed., 1999
[2] Barbero, "Introduction to Composite Materials Design", 3rd Ed., 2015

Based on FinEtoolsFlexStructures.jl by Petr Krysl - MIT License
=#

using LinearAlgebra: I

# =============================================================================
# Ply Definition
# =============================================================================

"""
    Ply(name, E1, E2, G12, ν12, thickness, angle; G13=G12, G23=G12, ρ=0.0)

Define a single ply (lamina) with orthotropic material properties.

# Arguments
- `name::String`: Ply identifier
- `E1`: Young's modulus in fiber direction [Pa]
- `E2`: Young's modulus transverse to fibers [Pa]  
- `G12`: In-plane shear modulus [Pa]
- `ν12`: Major Poisson's ratio [-]
- `thickness`: Ply thickness [m]
- `angle`: Fiber orientation angle from laminate x-axis [degrees]

# Optional
- `G13`: Transverse shear modulus in 1-3 plane (default: G12) [Pa]
- `G23`: Transverse shear modulus in 2-3 plane (default: G12) [Pa]
- `ρ`: Mass density [kg/m³]

# Example
```julia
ply = Ply("MyMaterial", 140e9, 10e9, 5e9, 0.3, 0.001, 45.0; ρ=1600.0)
```
"""
struct Ply
    name::String
    E1::Float64      # Fiber direction modulus [Pa]
    E2::Float64      # Transverse modulus [Pa]
    G12::Float64     # In-plane shear modulus [Pa]
    ν12::Float64     # Major Poisson's ratio
    G13::Float64     # Transverse shear 1-3 [Pa]
    G23::Float64     # Transverse shear 2-3 [Pa]
    thickness::Float64  # [m]
    angle::Float64   # Fiber angle [degrees]
    ρ::Float64       # Density [kg/m³]
    
    # Computed: plane stress stiffness in ply coordinates (3×3)
    Qbar::Matrix{Float64}
    # Computed: transverse shear stiffness in ply coordinates (2×2)
    Qts::Matrix{Float64}
    
    function Ply(name::String, E1::Real, E2::Real, G12::Real, ν12::Real,
                 thickness::Real, angle::Real;
                 G13::Real=G12, G23::Real=G12, ρ::Real=0.0)
        
        @assert 0.0 < ν12 < 1.0 "Poisson's ratio must be in (0, 1)"
        @assert thickness > 0 "Thickness must be positive"
        
        # Minor Poisson's ratio from symmetry: ν21 = ν12 × E2/E1
        ν21 = ν12 * E2 / E1
        
        # Plane stress stiffness matrix Q (in ply coordinates)
        denom = 1.0 - ν12 * ν21
        Q11 = E1 / denom
        Q22 = E2 / denom
        Q12 = ν12 * E2 / denom
        Q66 = G12
        
        Q = [Q11 Q12 0.0;
             Q12 Q22 0.0;
             0.0 0.0 Q66]
        
        Qts = [G13 0.0; 0.0 G23]
        
        new(name, Float64(E1), Float64(E2), Float64(G12), Float64(ν12),
            Float64(G13), Float64(G23), Float64(thickness), Float64(angle),
            Float64(ρ), Q, Qts)
    end
end

"""
    isotropic_ply(name, E, ν, thickness; ρ=0.0)

Create an isotropic ply (E1 = E2, standard shear modulus).
"""
function isotropic_ply(name::String, E::Real, ν::Real, thickness::Real; ρ::Real=0.0)
    G = E / (2 * (1 + ν))
    return Ply(name, E, E, G, ν, thickness, 0.0; G13=G, G23=G, ρ=ρ)
end

# =============================================================================
# Transformation Matrices (internal)
# =============================================================================

function _stress_transform_matrix(θ::Float64)
    m, n = cos(θ), sin(θ)
    return [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2]
end

function _strain_transform_matrix(θ::Float64)
    m, n = cos(θ), sin(θ)
    return [m^2 n^2 m*n; n^2 m^2 -m*n; -2*m*n 2*m*n m^2-n^2]
end

function _transform_Q_to_laminate(Q::Matrix{Float64}, θ::Float64)
    T = _stress_transform_matrix(θ)
    Tbar = _strain_transform_matrix(θ)
    return T \ (Q * Tbar)
end

function _transform_Qts_to_laminate(Qts::Matrix{Float64}, θ::Float64)
    m, n = cos(θ), sin(θ)
    T_ts = [m -n; n m]
    return T_ts' * Qts * T_ts
end

# =============================================================================
# Laminate Definition
# =============================================================================

"""
    Laminate(name, plies; offset=0.0)

Define a composite laminate from a stack of plies.

Plies are stacked from bottom (z = -t/2) to top (z = +t/2).

# Arguments
- `name::String`: Laminate identifier
- `plies::Vector{Ply}`: Stack of plies (first = bottom)
- `offset::Float64`: Offset of reference surface from mid-plane [m]

# Example
```julia
ply_0 = Ply("Mat", 100e9, 10e9, 5e9, 0.3, 0.001, 0.0)
ply_90 = Ply("Mat", 100e9, 10e9, 5e9, 0.3, 0.001, 90.0)
laminate = Laminate("MyLaminate", [ply_0, ply_90, ply_90, ply_0])
```
"""
struct Laminate
    name::String
    plies::Vector{Ply}
    offset::Float64
    total_thickness::Float64
    
    function Laminate(name::String, plies::Vector{Ply}; offset::Real=0.0)
        t_total = sum(p.thickness for p in plies)
        new(name, plies, Float64(offset), t_total)
    end
end

"""Total thickness of the laminate."""
thickness(lam::Laminate) = lam.total_thickness

# =============================================================================
# ABD Matrix Computation
# =============================================================================

"""
    laminate_stiffnesses(lam::Laminate) -> (A, B, D)

Compute laminate stiffness matrices using Classical Laminate Theory.

Returns:
- `A`: Membrane stiffness [N/m] (3×3)
- `B`: Extension-bending coupling [N] (3×3)
- `D`: Bending stiffness [N·m] (3×3)

The constitutive relation is:
```
[N]   [A  B] [ε₀]
[M] = [B  D] [κ ]
```
"""
function laminate_stiffnesses(lam::Laminate)
    A, B, D = zeros(3, 3), zeros(3, 3), zeros(3, 3)
    
    t_total = lam.total_thickness
    z = -t_total / 2 - lam.offset
    
    for ply in lam.plies
        z_bot, z_top = z, z + ply.thickness
        θ = ply.angle * π / 180
        Qbar = _transform_Q_to_laminate(ply.Qbar, θ)
        
        A .+= Qbar .* (z_top - z_bot)
        B .+= Qbar .* (z_top^2 - z_bot^2) / 2
        D .+= Qbar .* (z_top^3 - z_bot^3) / 3
        
        z = z_top
    end
    
    return A, B, D
end

"""
    laminate_transverse_shear_stiffness(lam::Laminate) -> H

Compute transverse shear stiffness with Vinson-Sierakowski correction.

Returns 2×2 matrix H for: [Qxz, Qyz]' = H × [γxz, γyz]'
"""
function laminate_transverse_shear_stiffness(lam::Laminate)
    H = zeros(2, 2)
    t_total = lam.total_thickness
    z = -t_total / 2 - lam.offset
    
    for ply in lam.plies
        z_bot, z_top = z, z + ply.thickness
        θ = ply.angle * π / 180
        Qts_bar = _transform_Qts_to_laminate(ply.Qts, θ)
        
        H .+= (5.0/4.0) .* Qts_bar .* (
            (z_top - z_bot) - (4.0/3.0) * (z_top^3 - z_bot^3) / t_total^2
        )
        z = z_top
    end
    
    return H
end

"""
    laminate_inertias(lam::Laminate) -> (ρ_area, ρ_inertia)

Compute mass per unit area [kg/m²] and rotational inertia density [kg].
"""
function laminate_inertias(lam::Laminate)
    ρ_area, ρ_inertia = 0.0, 0.0
    t_total = lam.total_thickness
    z = -t_total / 2 - lam.offset
    
    for ply in lam.plies
        z_bot, z_top = z, z + ply.thickness
        ρ_area += ply.ρ * (z_top - z_bot)
        ρ_inertia += ply.ρ * (z_top^3 - z_bot^3) / 3
        z = z_top
    end
    
    return ρ_area, ρ_inertia
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    symmetric_laminate(name, half_plies)

Create a symmetric laminate from half the ply stack.
Example: [0/45] → [0/45/45/0]
"""
function symmetric_laminate(name::String, half_plies::Vector{Ply})
    return Laminate(name, vcat(half_plies, reverse(half_plies)))
end
