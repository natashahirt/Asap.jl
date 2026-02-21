#=
Shell Elements for Asap.jl
==========================

Full shell element implementation with membrane, bending, and transverse shear stiffness.
Based on the T3FF (Triangular Flat-Facet) formulation with Discrete Shear Gap (DSG) technology.

Mathematical foundation:
- Membrane: constant strain triangle (CST)
- Bending: Discrete Kirchhoff Triangle (DKT) style curvature-displacement
- Transverse Shear: DSG formulation to avoid shear locking

References:
[1] Cui et al, "Analysis of plates and shells using an edge-based smoothed FEM"
    Comput Mech (2010) 45:141–156, DOI 10.1007/s00466-009-0429-9
[2] Lyly, Stenberg, Vihinen, "A stable bilinear element for Reissner-Mindlin plate"
    CMAME 110 (1993) 343-357

Original FinEtools implementation by Petr Krysl (FinEtoolsFlexStructures.jl)
Adapted for Asap.jl with permission - MIT License
=#

using LinearAlgebra: norm, dot, cross, mul!, lmul!

# Inline mean to avoid Statistics dependency
_mean(x) = sum(x) / length(x)

# =============================================================================
# Abstract Types
# =============================================================================

"""Abstract type for all shell elements."""
abstract type ShellElement <: AbstractElement end

# =============================================================================
# ShellTri3: 3-Node Triangular Shell Element
# =============================================================================

"""
    ShellTri3(nodes, thickness, E, ν; id=:shell, ρ=0.0, κ=5/6, drilling_scale=1.0, shear_stab=5/18)

A 3-node triangular shell element with full membrane, bending, and transverse shear stiffness.

Uses 6 DOFs per node: [u, v, w, θx, θy, θz] for a total of 18 DOFs.
Suitable for flat or nearly-flat shell structures, plates, and slabs.

# Arguments
- `nodes::NTuple{3, Node}`: Three corner nodes (counter-clockwise order defines +Z normal)
- `thickness::Quantity{Length}`: Shell thickness
- `E::Quantity{Pressure}`: Young's modulus
- `ν::Real`: Poisson's ratio

# Optional
- `id::Symbol`: Element identifier (default: `:shell`)
- `ρ::Real`: Mass density in kg/m³ (default: 0.0)
- `κ::Real`: Shear correction factor (default: 5/6 ≈ 0.833 for rectangular sections)
- `drilling_scale::Real`: Drilling DOF stabilization scale (default: 1.0)
- `shear_stab::Real`: Lyly-Stenberg-Vihinen shear stabilization factor (default: 5/18 ≈ 0.278)

# Formulation Parameters

**Shear correction factor (κ)**: Accounts for non-uniform shear stress distribution through 
thickness. Standard value 5/6 is for rectangular cross-sections. Use:
- κ = 5/6 (≈0.833): Rectangular sections (default)
- κ = 0.9: Parabolic shear distribution approximation
- κ = 1.0: No correction (conservative for thin shells)

**Drilling scale**: Scales the artificial drilling DOF (θz) stiffness relative to average
bending stiffness. Values 0.1–1.0 are typical. Set to 0.0 to disable (not recommended).

**Shear stabilization (shear_stab)**: Controls the Lyly-Stenberg-Vihinen transverse shear 
stabilization. Factor appears as: t³/(t² + shear_stab·h²) where h² ≈ 2·Area.
Standard value 5/18 from [Lyly, Stenberg, Vihinen 1993].

# Example
```julia
n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
n3 = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)

# Default parameters
shell = ShellTri3((n1, n2, n3), 0.2u"m", 30e9u"Pa", 0.2)

# Custom parameters for thick shell
shell = ShellTri3((n1, n2, n3), 0.5u"m", 30e9u"Pa", 0.2; κ=0.9, drilling_scale=0.5)
```

# Mathematical Formulation
The element stiffness is K = Kₘ + Kᵦ + Kₛ + Kᵈ where:
- Kₘ: Membrane stiffness (in-plane)
- Kᵦ: Bending stiffness (curvature)
- Kₛ: Transverse shear stiffness (DSG formulation)
- Kᵈ: Drilling DOF stabilization

The shear stabilization factor follows Lyly-Stenberg-Vihinen (1993).
"""
mutable struct ShellTri3 <: ShellElement
    nodes::NTuple{3, Node}
    thickness::Float64          # [m]
    E::Float64                  # [Pa]
    ν::Float64                  # [-]
    ρ::Float64                  # [kg/m³]
    κ::Float64                  # Shear correction factor
    drilling_scale::Float64     # Drilling DOF stabilization scale
    shear_stab::Float64         # Shear stabilization factor
    elementID::Int64
    globalID::Vector{Int64}     # 18 DOFs
    area::Float64               # [m²]
    K::Matrix{Float64}          # 18×18 global stiffness
    M::Matrix{Float64}          # 18×18 lumped mass (diagonal)
    R::Matrix{Float64}          # 18×18 transformation GCS→LCS
    LCS::Vector{Vector{Float64}}  # Local axes [x̂, ŷ, ẑ] in GCS
    ecoords_e::Matrix{Float64}  # Node coords in element frame (3×2)
    id::Symbol
    
    function ShellTri3(
        nodes::NTuple{3, Node},
        thickness::Quantity,
        E::Quantity,
        ν::Real;
        id::Symbol = :shell,
        ρ::Real = 0.0,
        κ::Real = 5.0/6.0,
        drilling_scale::Real = 1.0,
        shear_stab::Real = 5.0/18.0
    )
        @assert 0.0 <= ν < 0.5 "Poisson's ratio must be in [0, 0.5)"
        @assert 0.0 < κ <= 1.0 "Shear correction factor must be in (0, 1]"
        @assert drilling_scale >= 0.0 "Drilling scale must be non-negative"
        @assert shear_stab >= 0.0 "Shear stabilization must be non-negative"
        
        t = ustrip(u"m", thickness)
        E_val = ustrip(u"Pa", E)
        area = _compute_triangle_area(nodes)
        
        new(
            nodes,
            t,
            E_val,
            Float64(ν),
            Float64(ρ),
            Float64(κ),
            Float64(drilling_scale),
            Float64(shear_stab),
            0,
            Vector{Int64}(undef, 18),
            area,
            zeros(18, 18),
            zeros(18, 18),
            zeros(18, 18),
            [zeros(3), zeros(3), zeros(3)],
            zeros(3, 2),
            id
        )
    end
    
    # Constructor with ShellMaterial
    function ShellTri3(
        nodes::NTuple{3, Node},
        thickness::Quantity,
        material::ShellMaterial;
        id::Symbol = :shell,
        κ::Real = 5.0/6.0,
        drilling_scale::Real = 1.0,
        shear_stab::Real = 5.0/18.0
    )
        @assert 0.0 < κ <= 1.0 "Shear correction factor must be in (0, 1]"
        @assert drilling_scale >= 0.0 "Drilling scale must be non-negative"
        @assert shear_stab >= 0.0 "Shear stabilization must be non-negative"
        
        t = ustrip(u"m", thickness)
        area = _compute_triangle_area(nodes)
        
        new(
            nodes,
            t,
            material.E,
            material.ν,
            material.ρ,
            Float64(κ),
            Float64(drilling_scale),
            Float64(shear_stab),
            0,
            Vector{Int64}(undef, 18),
            area,
            zeros(18, 18),
            zeros(18, 18),
            zeros(18, 18),
            [zeros(3), zeros(3), zeros(3)],
            zeros(3, 2),
            id
        )
    end
    
    # Constructor with Material (Asap's frame material type)
    function ShellTri3(
        nodes::NTuple{3, Node},
        thickness::Quantity,
        material::Material;
        id::Symbol = :shell
    )
        t = ustrip(u"m", thickness)
        E_val = ustrip(u"Pa", material.E)
        ρ_val = ustrip(u"kg/m^3", material.ρ)
        area = _compute_triangle_area(nodes)
        
        new(
            nodes,
            t,
            E_val,
            material.ν,
            ρ_val,
            0,
            Vector{Int64}(undef, 18),
            area,
            zeros(18, 18),
            zeros(18, 18),
            zeros(18, 18),
            [zeros(3), zeros(3), zeros(3)],
            zeros(3, 2),
            id
        )
    end
end

# =============================================================================
# Core Geometry Functions
# =============================================================================

"""Extract node coordinates as 3×3 matrix (row = node, col = x,y,z)."""
function _get_coords(elem::ShellTri3)
    ecoords = zeros(3, 3)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:3
            ecoords[i, j] = ustrip(u"m", node.position[j])
        end
    end
    return ecoords
end

"""Compute triangle area from three nodes: Area = 0.5 * ||(p2-p1) × (p3-p1)||"""
function _compute_triangle_area(nodes::NTuple{3, Node})
    p1 = [ustrip(u"m", nodes[1].position[j]) for j in 1:3]
    p2 = [ustrip(u"m", nodes[2].position[j]) for j in 1:3]
    p3 = [ustrip(u"m", nodes[3].position[j]) for j in 1:3]
    v1 = p2 - p1
    v2 = p3 - p1
    cross_prod = [v1[2]*v2[3] - v1[3]*v2[2],
                  v1[3]*v2[1] - v1[1]*v2[3],
                  v1[1]*v2[2] - v1[2]*v2[1]]
    return 0.5 * sqrt(sum(c^2 for c in cross_prod))
end

"""Compute Jacobian matrix J0 (3×2) from node coordinates."""
function _compute_J!(J0::Matrix{Float64}, ecoords::Matrix{Float64})
    # J0[:,1] = edge 1→2
    # J0[:,2] = edge 1→3
    for j in 1:3
        J0[j, 1] = ecoords[2, j] - ecoords[1, j]
        J0[j, 2] = ecoords[3, j] - ecoords[1, j]
    end
    return J0
end

"""
Compute element local coordinate system E_G (3×3).

Columns are: [x̂_e, ŷ_e, ẑ_e] expressed in global coordinates.
- x̂_e: along edge 1→2 (normalized)
- ẑ_e: normal to element (J0[:,1] × J0[:,2], normalized)
- ŷ_e: ẑ_e × x̂_e (completes right-handed system)
"""
function _compute_E_G!(E_G::Matrix{Float64}, J0::Matrix{Float64})
    # x̂ = J0[:,1] normalized
    for j in 1:3
        E_G[j, 1] = J0[j, 1]
    end
    n = sqrt(E_G[1,1]^2 + E_G[2,1]^2 + E_G[3,1]^2)
    for j in 1:3
        E_G[j, 1] /= n
    end
    
    # ẑ = x̂ × J0[:,2] normalized (using J0[:,1] not normalized x̂)
    E_G[1, 3] = -E_G[3, 1] * J0[2, 2] + E_G[2, 1] * J0[3, 2]
    E_G[2, 3] =  E_G[3, 1] * J0[1, 2] - E_G[1, 1] * J0[3, 2]
    E_G[3, 3] = -E_G[2, 1] * J0[1, 2] + E_G[1, 1] * J0[2, 2]
    n = sqrt(E_G[1,3]^2 + E_G[2,3]^2 + E_G[3,3]^2)
    for j in 1:3
        E_G[j, 3] /= n
    end
    
    # ŷ = ẑ × x̂
    E_G[1, 2] = -E_G[3, 3] * E_G[2, 1] + E_G[2, 3] * E_G[3, 1]
    E_G[2, 2] =  E_G[3, 3] * E_G[1, 1] - E_G[1, 3] * E_G[3, 1]
    E_G[3, 2] = -E_G[2, 3] * E_G[1, 1] + E_G[1, 3] * E_G[2, 1]
    
    return E_G
end

"""Project global coords to element local 2D coords (3×2 matrix)."""
function _compute_ecoords_e!(ecoords_e::Matrix{Float64}, J0::Matrix{Float64}, E_G::Matrix{Float64})
    # Node 1 is at origin
    ecoords_e[1, 1] = 0.0
    ecoords_e[1, 2] = 0.0
    # Node 2: project J0[:,1] onto x̂, ŷ
    ecoords_e[2, 1] = dot(view(J0, :, 1), view(E_G, :, 1))
    ecoords_e[2, 2] = dot(view(J0, :, 1), view(E_G, :, 2))
    # Node 3: project J0[:,2] onto x̂, ŷ  
    ecoords_e[3, 1] = dot(view(J0, :, 2), view(E_G, :, 1))
    ecoords_e[3, 2] = dot(view(J0, :, 2), view(E_G, :, 2))
    return ecoords_e
end

"""Compute shape function gradients and element area."""
function _compute_gradN_Ae!(gradN_e::Matrix{Float64}, ecoords_e::Matrix{Float64})
    # Triangle with vertices at (x1,y1), (x2,y2), (x3,y3) in local coords
    # Using standard area coordinates
    a = ecoords_e[2, 1] - ecoords_e[1, 1]  # x2 - x1
    b = ecoords_e[2, 2] - ecoords_e[1, 2]  # y2 - y1
    c = ecoords_e[3, 1] - ecoords_e[1, 1]  # x3 - x1
    d = ecoords_e[3, 2] - ecoords_e[1, 2]  # y3 - y1
    
    J = a * d - b * c  # 2 × Area
    
    # ∂N/∂x, ∂N/∂y for each node
    gradN_e[1, 1] = (b - d) / J
    gradN_e[2, 1] = d / J
    gradN_e[3, 1] = -b / J
    gradN_e[1, 2] = (c - a) / J
    gradN_e[2, 2] = -c / J
    gradN_e[3, 2] = a / J
    
    return gradN_e, J / 2.0  # Area
end

# =============================================================================
# Strain-Displacement Matrices (B-matrices)
# =============================================================================

"""
Membrane strain-displacement matrix Bm (3×18).

Maps nodal DOFs to in-plane strains [εxx, εyy, γxy].
Uses constant strain triangle (CST) formulation.
"""
function _Bmmat!(Bm::Matrix{Float64}, gradN::Matrix{Float64})
    fill!(Bm, 0.0)
    for i in 1:3
        offset = (i - 1) * 6
        Bm[1, offset + 1] = gradN[i, 1]       # ∂u/∂x → εxx
        Bm[2, offset + 2] = gradN[i, 2]       # ∂v/∂y → εyy
        Bm[3, offset + 1] = gradN[i, 2]       # ∂u/∂y → γxy
        Bm[3, offset + 2] = gradN[i, 1]       # ∂v/∂x → γxy
    end
    return Bm
end

"""
Bending curvature-displacement matrix Bb (3×18).

Maps nodal DOFs to curvatures [κxx, κyy, κxy].
Uses DKT-style formulation.
"""
function _Bbmat!(Bb::Matrix{Float64}, gradN::Matrix{Float64})
    fill!(Bb, 0.0)
    for i in 1:3
        offset = (i - 1) * 6
        Bb[1, offset + 5] = gradN[i, 1]       # ∂θy/∂x → κxx
        Bb[2, offset + 4] = -gradN[i, 2]      # -∂θx/∂y → κyy
        Bb[3, offset + 4] = -gradN[i, 1]      # -∂θx/∂x → κxy
        Bb[3, offset + 5] = gradN[i, 2]       # ∂θy/∂y → κxy
    end
    return Bb
end

"""
Transverse shear strain-displacement matrix Bs (2×18).

Uses Discrete Shear Gap (DSG) formulation to avoid shear locking.
Averages contributions from three node orderings.
"""
function _Bsmat!(Bs::Matrix{Float64}, ecoords_e::Matrix{Float64}, Ae::Float64)
    fill!(Bs, 0.0)
    
    # Average three orderings for isotropy
    for ordering in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
        _add_Bsmat_ordering!(Bs, ecoords_e, Ae, ordering)
    end
    
    Bs .*= (1.0 / 3.0)
    return Bs
end

"""Add DSG shear contribution for one node ordering."""
function _add_Bsmat_ordering!(Bs::Matrix{Float64}, ecoords_e::Matrix{Float64}, 
                              Ae::Float64, ordering::NTuple{3, Int})
    s, p, q = ordering
    
    a = ecoords_e[p, 1] - ecoords_e[s, 1]
    b = ecoords_e[p, 2] - ecoords_e[s, 2]
    c = ecoords_e[q, 1] - ecoords_e[s, 1]
    d = ecoords_e[q, 2] - ecoords_e[s, 2]
    
    m = 1.0 / (2.0 * Ae)
    
    # Node s contributions
    co = (s - 1) * 6
    Bs[1, co + 3] += m * (b - d)
    Bs[1, co + 5] += m * Ae
    Bs[2, co + 3] += m * (c - a)
    Bs[2, co + 4] += m * (-Ae)
    
    # Node p contributions
    co = (p - 1) * 6
    Bs[1, co + 3] += m * d
    Bs[1, co + 4] += m * (-b * d / 2)
    Bs[1, co + 5] += m * (a * d / 2)
    Bs[2, co + 3] += m * (-c)
    Bs[2, co + 4] += m * (b * c / 2)
    Bs[2, co + 5] += m * (-a * c / 2)
    
    # Node q contributions
    co = (q - 1) * 6
    Bs[1, co + 3] += m * (-b)
    Bs[1, co + 4] += m * (b * d / 2)
    Bs[1, co + 5] += m * (-b * c / 2)
    Bs[2, co + 3] += m * a
    Bs[2, co + 4] += m * (-a * d / 2)
    Bs[2, co + 5] += m * (a * c / 2)
    
    return Bs
end

# =============================================================================
# Material Stiffness Matrices
# =============================================================================

"""
Plane stress constitutive matrix Dps (3×3).

For membrane and bending: σ = Dps · ε
"""
function _plane_stress_D(E::Float64, ν::Float64)
    c = E / (1.0 - ν^2)
    return [c      c*ν    0.0;
            c*ν    c      0.0;
            0.0    0.0    c*(1.0-ν)/2.0]
end

"""
Transverse shear constitutive matrix Dt (2×2).

For transverse shear: τ = Dt · γ
"""
function _transverse_shear_D(E::Float64, ν::Float64, κ::Float64)
    G = E / (2.0 * (1.0 + ν))
    return [κ*G  0.0;
            0.0  κ*G]
end

# =============================================================================
# Stiffness Matrix Assembly
# =============================================================================

"""Add B'DB contribution to upper triangle of K, then complete lower triangle."""
function _add_btdb!(K::Matrix{Float64}, B::Matrix{Float64}, 
                    factor::Float64, D::Matrix{Float64}, DB::Matrix{Float64})
    # DB = D * B
    mul!(DB, D, B)
    
    # K += factor * B' * DB (upper triangle only)
    n = size(B, 2)
    for j in 1:n
        for i in 1:j
            for k in axes(B, 1)
                K[i, j] += factor * B[k, i] * DB[k, j]
            end
        end
    end
    
    return K
end

"""Complete lower triangle from upper triangle (symmetric matrix)."""
function _complete_lt!(K::Matrix{Float64})
    n = size(K, 1)
    for i in 2:n
        for j in 1:(i-1)
            K[i, j] = K[j, i]
        end
    end
    return K
end

"""
    local_K(elem::ShellTri3) -> Matrix{Float64}

Compute 18×18 element stiffness matrix in local coordinates.

The stiffness combines:
- Membrane (in-plane stretch/shear)
- Bending (plate curvature)
- Transverse shear (DSG formulation with Lyly-Stenberg-Vihinen stabilization)
- Drilling DOF stabilization

Uses element's configurable parameters: κ (shear correction), drilling_scale, shear_stab.
"""
function local_K(elem::ShellTri3)
    t = elem.thickness
    E = elem.E
    ν = elem.ν
    κ = elem.κ
    Ae = elem.area
    ecoords_e = elem.ecoords_e
    
    # Preallocate
    K = zeros(18, 18)
    gradN_e = zeros(3, 2)
    Bm = zeros(3, 18)
    Bb = zeros(3, 18)
    Bs = zeros(2, 18)
    DBm = zeros(3, 18)
    DBb = zeros(3, 18)
    DBs = zeros(2, 18)
    
    # Compute gradients
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, ecoords_e)
    
    # Material matrices
    Dps = _plane_stress_D(E, ν)
    Dt = _transverse_shear_D(E, ν, κ)
    
    # Membrane stiffness: K_m = t * Ae * Bm' * Dps * Bm
    _Bmmat!(Bm, gradN_e)
    _add_btdb!(K, Bm, t * Ae, Dps, DBm)
    
    # Bending stiffness: K_b = (t³/12) * Ae * Bb' * Dps * Bb
    _Bbmat!(Bb, gradN_e)
    _add_btdb!(K, Bb, (t^3 / 12.0) * Ae, Dps, DBb)
    
    # Transverse shear stiffness with Lyly-Stenberg-Vihinen stabilization
    # Shear factor: t³/(t² + shear_stab·h²) where h² ≈ 2·Ae
    shear_factor = (t^3 / (t^2 + elem.shear_stab * 2.0 * Ae)) * Ae
    
    _Bsmat!(Bs, ecoords_e, Ae)
    _add_btdb!(K, Bs, shear_factor, Dt, DBs)
    
    # Complete symmetric matrix
    _complete_lt!(K)
    
    # Drilling DOF stabilization (θz)
    # Add stiffness to drilling DOFs based on average bending stiffness
    if elem.drilling_scale > 0.0
        kavg = _mean([K[4,4], K[5,5], K[10,10], K[11,11], K[16,16], K[17,17]])
        K[6, 6] += kavg * elem.drilling_scale
        K[12, 12] += kavg * elem.drilling_scale
        K[18, 18] += kavg * elem.drilling_scale
    end
    
    return K
end

# =============================================================================
# Transformation Matrix
# =============================================================================

"""
Build 18×18 transformation matrix from global to local coordinates.

For each node, transforms [u,v,w,θx,θy,θz]_global → [u,v,w,θx,θy,θz]_local.
"""
function _build_transformation!(R::Matrix{Float64}, E_G::Matrix{Float64})
    fill!(R, 0.0)
    
    # E_G columns are local basis vectors in global coords
    # To transform from global to local: multiply by E_G'
    T3 = E_G'  # 3×3 transformation
    
    for node_idx in 1:3
        offset = (node_idx - 1) * 6
        
        # Translation DOFs
        for i in 1:3
            for j in 1:3
                R[offset + i, offset + j] = T3[i, j]
            end
        end
        
        # Rotation DOFs (same transformation)
        for i in 1:3
            for j in 1:3
                R[offset + 3 + i, offset + 3 + j] = T3[i, j]
            end
        end
    end
    
    return R
end

# =============================================================================
# Element Processing Functions
# =============================================================================

"""
    lcs!(elem::ShellTri3)

Compute local coordinate system and element geometry.
"""
function lcs!(elem::ShellTri3)
    ecoords = _get_coords(elem)
    
    J0 = zeros(3, 2)
    E_G = zeros(3, 3)
    
    _compute_J!(J0, ecoords)
    _compute_E_G!(E_G, J0)
    _compute_ecoords_e!(elem.ecoords_e, J0, E_G)
    
    # Store LCS axes
    elem.LCS[1] = E_G[:, 1]
    elem.LCS[2] = E_G[:, 2]
    elem.LCS[3] = E_G[:, 3]
    
    # Compute area
    gradN = zeros(3, 2)
    _, elem.area = _compute_gradN_Ae!(gradN, elem.ecoords_e)
    
    return elem
end

"""
    R!(elem::ShellTri3)

Compute and store transformation matrix.
"""
function R!(elem::ShellTri3)
    E_G = hcat(elem.LCS...)  # Reconstruct E_G from stored LCS
    _build_transformation!(elem.R, E_G)
    return elem
end

"""
    global_K!(elem::ShellTri3)

Compute and store global stiffness matrix.
"""
function global_K!(elem::ShellTri3)
    K_local = local_K(elem)
    elem.K = elem.R' * K_local * elem.R
    return elem
end

"""
    process!(elem::ShellTri3)

Full element processing: compute LCS, R, K, and M.
"""
function process!(elem::ShellTri3)
    lcs!(elem)
    R!(elem)
    global_K!(elem)
    
    # Lumped mass matrix (if density specified)
    if elem.ρ > 0.0
        mass_per_node = elem.ρ * elem.thickness * elem.area / 3.0
        rot_inertia = elem.ρ * (elem.thickness^3 / 12.0) * elem.area / 3.0
        
        fill!(elem.M, 0.0)
        for i in 1:3
            offset = (i - 1) * 6
            # Translational mass
            for d in 1:3
                elem.M[offset + d, offset + d] = mass_per_node
            end
            # Rotational inertia
            for d in 4:6
                elem.M[offset + d, offset + d] = rot_inertia
            end
        end
    end
    
    return elem
end

# =============================================================================
# Stress Recovery
# =============================================================================

"""
    stress(elem::ShellTri3, u_global::Vector{Float64}) -> Vector{Float64}

Compute membrane stresses [σxx, σyy, τxy] in local coordinates.

# Arguments
- `elem`: Processed shell element
- `u_global`: Global displacement vector (must contain all 18 DOFs for this element)
"""
function stress(elem::ShellTri3, u_global::Vector{Float64})
    # Extract element DOFs from global vector
    u_elem = zeros(18)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:3  # Only translation DOFs for membrane
            global_dof = node.globalID[j]
            u_elem[(i-1)*6 + j] = u_global[global_dof]
        end
    end
    
    # Transform to local
    u_local = elem.R * u_elem
    
    # Compute membrane strain
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    Bm = zeros(3, 18)
    _Bmmat!(Bm, gradN_e)
    
    ε = Bm * u_local
    
    # Compute stress
    Dps = _plane_stress_D(elem.E, elem.ν)
    σ = Dps * ε
    
    return σ
end

"""
    bending_moments(elem::ShellTri3, u_global::Vector{Float64}) -> Vector{Float64}

Compute bending moments [Mxx, Myy, Mxy] in local coordinates.
"""
function bending_moments(elem::ShellTri3, u_global::Vector{Float64})
    # Extract element DOFs
    u_elem = zeros(18)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:6
            global_dof = node.globalID[j]
            u_elem[(i-1)*6 + j] = u_global[global_dof]
        end
    end
    
    # Transform to local
    u_local = elem.R * u_elem
    
    # Compute curvature
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    Bb = zeros(3, 18)
    _Bbmat!(Bb, gradN_e)
    
    κ = Bb * u_local
    
    # Compute moments: M = D_bending * κ where D_bending = (t³/12) * Dps
    Dps = _plane_stress_D(elem.E, elem.ν)
    M = (elem.thickness^3 / 12.0) * Dps * κ
    
    return M
end

# =============================================================================
# Zero-Allocation Bending Moments Workspace
# =============================================================================

"""
Thread-local workspace for zero-allocation `bending_moments!` calls.
Holds all temporary buffers needed to compute shell bending moments.
"""
struct ShellMomentWorkspace
    u_elem::Vector{Float64}     # length 18
    u_local::Vector{Float64}    # length 18
    gradN_e::Matrix{Float64}    # 3×2
    Bb::Matrix{Float64}         # 3×18
    κ::Vector{Float64}          # length 3
    Dps::Matrix{Float64}        # 3×3
end

function ShellMomentWorkspace()
    ShellMomentWorkspace(
        zeros(18), zeros(18), zeros(3, 2),
        zeros(3, 18), zeros(3), zeros(3, 3)
    )
end

const _MOMENT_WS = ShellMomentWorkspace[ShellMomentWorkspace() for _ in 1:Threads.nthreads()]

"""Fill plane-stress D matrix in-place."""
function _plane_stress_D!(D::Matrix{Float64}, E::Float64, ν::Float64)
    c = E / (1.0 - ν^2)
    D[1,1] = c;       D[1,2] = c*ν;     D[1,3] = 0.0
    D[2,1] = c*ν;     D[2,2] = c;       D[2,3] = 0.0
    D[3,1] = 0.0;     D[3,2] = 0.0;     D[3,3] = c*(1.0-ν)/2.0
    return D
end

"""
    bending_moments!(M, elem, u_global [, ws])

Compute bending moments `[Mxx, Myy, Mxy]` in local coordinates,
writing into the pre-allocated 3-vector `M`.  **Zero heap allocations.**
"""
function bending_moments!(M::AbstractVector{Float64}, elem::ShellTri3,
                          u_global::Vector{Float64},
                          ws::ShellMomentWorkspace = _MOMENT_WS[Threads.threadid()])
    # Extract element DOFs
    fill!(ws.u_elem, 0.0)
    @inbounds for i in 1:3
        base = (i - 1) * 6
        node = elem.nodes[i]
        for j in 1:6
            ws.u_elem[base + j] = u_global[node.globalID[j]]
        end
    end

    # Transform to local: u_local = R * u_elem
    mul!(ws.u_local, elem.R, ws.u_elem)

    # Curvature: κ = Bb * u_local
    _compute_gradN_Ae!(ws.gradN_e, elem.ecoords_e)
    _Bbmat!(ws.Bb, ws.gradN_e)
    mul!(ws.κ, ws.Bb, ws.u_local)

    # M = (t³/12) * Dps * κ
    _plane_stress_D!(ws.Dps, elem.E, elem.ν)
    mul!(M, ws.Dps, ws.κ)
    lmul!(elem.thickness^3 / 12.0, M)
    return M
end

# =============================================================================
# DOF Indexing for Preprocessing
# =============================================================================

"""Assign global DOF indices for a ShellTri3 element."""
function populate_globalID!(elem::ShellTri3)
    elem.globalID = vcat([n.globalID for n in elem.nodes]...)
    return elem
end