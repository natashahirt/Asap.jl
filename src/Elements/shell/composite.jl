#=
Composite Shell Elements for Asap.jl
====================================

Triangular shell element for layered (composite) materials using Classical 
Laminate Theory (CLT) and ABD matrix formulation.

Uses the same T3FF (Triangular Flat-Facet) geometry as ShellTri3 but with
laminate-based material stiffness instead of isotropic.

Mathematical foundation:
- ABD matrices computed from ply stack using CLT
- Membrane-bending coupling (B matrix) included
- Transverse shear with Vinson-Sierakowski correction

References:
[1] R.M. Jones, "Mechanics of Composite Materials", 2nd Ed., 1999
[2] Barbero, "Introduction to Composite Materials Design", 3rd Ed., 2015
[3] Cui et al, "Analysis of plates and shells using ES-FEM", Comput Mech 2010

Original FinEtools implementation by Petr Krysl (FinEtoolsFlexStructures.jl)
Adapted for Asap.jl with permission - MIT License
=#

# _mean is defined in shell.jl and shared across the module

# =============================================================================
# CompositeShellTri3: 3-Node Triangular Shell with Laminate Material
# =============================================================================

"""
    CompositeShellTri3(nodes, laminate; id=:composite_shell)

A 3-node triangular shell element for layered composite materials.

Uses 6 DOFs per node: [u, v, w, θx, θy, θz] for a total of 18 DOFs.
Material stiffness is computed from Classical Laminate Theory using the 
ABD matrix formulation.

# Arguments
- `nodes::NTuple{3, Node}`: Three corner nodes (counter-clockwise order defines +Z normal)
- `laminate::Laminate`: Composite laminate definition (see `Laminate`)

# Optional
- `id::Symbol`: Element identifier (default: `:composite_shell`)

# Example
```julia
# Define plies
ply_0 = Ply("T300", 181e9, 10.3e9, 7.17e9, 0.28, 0.125e-3, 0.0; ρ=1600.0)
ply_90 = Ply("T300", 181e9, 10.3e9, 7.17e9, 0.28, 0.125e-3, 90.0; ρ=1600.0)

# Create symmetric laminate [0/90/90/0]
laminate = Laminate("CFRP_sym", [ply_0, ply_90, ply_90, ply_0])

# Create nodes
n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
n3 = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)

# Create composite shell element
shell = CompositeShellTri3((n1, n2, n3), laminate)
```

# Constitutive Relation
The laminate constitutive relation is:
```
[N]   [A  B] [ε₀]
[M] = [B  D] [κ ]
```
where:
- N: membrane forces [Nxx, Nyy, Nxy] (N/m)
- M: bending moments [Mxx, Myy, Mxy] (N·m/m)
- ε₀: mid-plane strains
- κ: curvatures
- A: membrane stiffness (N/m)
- B: extension-bending coupling (N)
- D: bending stiffness (N·m)
"""
mutable struct CompositeShellTri3 <: ShellElement
    nodes::NTuple{3, Node}
    laminate::Laminate
    elementID::Int64
    globalID::Vector{Int64}     # 18 DOFs
    area::Float64               # [m²]
    K::Matrix{Float64}          # 18×18 global stiffness
    M::Matrix{Float64}          # 18×18 lumped mass
    R::Matrix{Float64}          # 18×18 transformation GCS→LCS
    LCS::Vector{Vector{Float64}}  # Local axes [x̂, ŷ, ẑ] in GCS
    ecoords_e::Matrix{Float64}  # Node coords in element frame (3×2)
    id::Symbol
    
    function CompositeShellTri3(
        nodes::NTuple{3, Node},
        laminate::Laminate;
        id::Symbol = :composite_shell
    )
        new(
            nodes,
            laminate,
            0,
            Vector{Int64}(undef, 18),
            0.0,
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
# Geometry Functions (shared with ShellTri3)
# =============================================================================

"""Extract node coordinates as 3×3 matrix."""
function _get_coords(elem::CompositeShellTri3)
    ecoords = zeros(3, 3)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:3
            ecoords[i, j] = ustrip(u"m", node.position[j])
        end
    end
    return ecoords
end

# =============================================================================
# Composite Stiffness Matrix Assembly
# =============================================================================

"""
    local_K(elem::CompositeShellTri3) -> Matrix{Float64}

Compute 18×18 element stiffness matrix in local coordinates.

Uses ABD matrices from Classical Laminate Theory:
- K = Bm'*A*Bm + Bm'*B*Bb + Bb'*B'*Bm + Bb'*D*Bb + Bs'*H*Bs + K_drilling

The extension-bending coupling (B matrix) allows proper modeling of 
asymmetric laminates.
"""
function local_K(elem::CompositeShellTri3)
    lam = elem.laminate
    Ae = elem.area
    ecoords_e = elem.ecoords_e
    t = thickness(lam)
    
    # Get ABD matrices
    A, B, D = laminate_stiffnesses(lam)
    H = laminate_transverse_shear_stiffness(lam)
    
    # Preallocate
    K = zeros(18, 18)
    gradN_e = zeros(3, 2)
    Bm = zeros(3, 18)
    Bb = zeros(3, 18)
    Bs = zeros(2, 18)
    
    # Temporary storage for matrix products
    ABm = zeros(3, 18)
    BBm = zeros(3, 18)
    BBb = zeros(3, 18)
    DBb = zeros(3, 18)
    HBs = zeros(2, 18)
    
    # Compute gradients
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, ecoords_e)
    
    # Build strain-displacement matrices
    _Bmmat!(Bm, gradN_e)
    _Bbmat!(Bb, gradN_e)
    _Bsmat!(Bs, ecoords_e, Ae)
    
    # --- Membrane stiffness: Bm' * A * Bm ---
    mul!(ABm, A, Bm)
    _add_btc!(K, Bm, ABm, Ae)
    
    # --- Extension-bending coupling: Bm' * B * Bb + Bb' * B' * Bm ---
    mul!(BBb, B, Bb)
    _add_btc!(K, Bm, BBb, Ae)  # Bm' * B * Bb
    
    mul!(BBm, B', Bm)  # B' * Bm (B is generally not symmetric)
    _add_btc!(K, Bb, BBm, Ae)  # Bb' * B' * Bm
    
    # --- Bending stiffness: Bb' * D * Bb ---
    mul!(DBb, D, Bb)
    _add_btc!(K, Bb, DBb, Ae)
    
    # --- Transverse shear stiffness with stabilization ---
    # Shear factor from Lyly-Stenberg-Vihinen stabilization
    mult_el_size = 5.0 / 12.0 / 1.5
    shear_factor = t^2 / (t^2 + mult_el_size * 2.0 * Ae)
    
    mul!(HBs, H, Bs)
    _add_btc!(K, Bs, HBs, shear_factor * Ae)
    
    # Make symmetric (numerical cleanup)
    for i in 2:18
        for j in 1:(i-1)
            avg = (K[i,j] + K[j,i]) / 2
            K[i,j] = avg
            K[j,i] = avg
        end
    end
    
    # --- Drilling DOF stabilization ---
    kavg = _mean([K[4,4], K[5,5], K[10,10], K[11,11], K[16,16], K[17,17]])
    drilling_scale = 1.0
    K[6, 6] += kavg * drilling_scale
    K[12, 12] += kavg * drilling_scale
    K[18, 18] += kavg * drilling_scale
    
    return K
end

"""Add B' * C contribution to K (where C = D*B for some stiffness D)."""
function _add_btc!(K::Matrix{Float64}, B::Matrix{Float64}, C::Matrix{Float64}, factor::Float64)
    # K += factor * B' * C
    n = size(B, 2)
    for j in 1:n
        for i in 1:n
            for k in axes(B, 1)
                K[i, j] += factor * B[k, i] * C[k, j]
            end
        end
    end
    return K
end

# =============================================================================
# Element Processing
# =============================================================================

"""
    lcs!(elem::CompositeShellTri3)

Compute local coordinate system and element geometry.
"""
function lcs!(elem::CompositeShellTri3)
    ecoords = _get_coords(elem)
    
    J0 = zeros(3, 2)
    E_G = zeros(3, 3)
    
    _compute_J!(J0, ecoords)
    _compute_E_G!(E_G, J0)
    _compute_ecoords_e!(elem.ecoords_e, J0, E_G)
    
    elem.LCS[1] = E_G[:, 1]
    elem.LCS[2] = E_G[:, 2]
    elem.LCS[3] = E_G[:, 3]
    
    gradN = zeros(3, 2)
    _, elem.area = _compute_gradN_Ae!(gradN, elem.ecoords_e)
    
    return elem
end

"""
    R!(elem::CompositeShellTri3)

Compute and store transformation matrix.
"""
function R!(elem::CompositeShellTri3)
    E_G = hcat(elem.LCS...)
    _build_transformation!(elem.R, E_G)
    return elem
end

"""
    global_K!(elem::CompositeShellTri3)

Compute and store global stiffness matrix.
"""
function global_K!(elem::CompositeShellTri3)
    K_local = local_K(elem)
    elem.K = elem.R' * K_local * elem.R
    return elem
end

"""
    process!(elem::CompositeShellTri3)

Full element processing: compute LCS, R, K, and M.
"""
function process!(elem::CompositeShellTri3)
    lcs!(elem)
    R!(elem)
    global_K!(elem)
    
    # Lumped mass matrix from laminate inertias
    ρ_area, ρ_inertia = laminate_inertias(elem.laminate)
    
    if ρ_area > 0.0
        mass_per_node = ρ_area * elem.area / 3.0
        rot_per_node = ρ_inertia * elem.area / 3.0
        
        fill!(elem.M, 0.0)
        for i in 1:3
            offset = (i - 1) * 6
            for d in 1:3
                elem.M[offset + d, offset + d] = mass_per_node
            end
            for d in 4:6
                elem.M[offset + d, offset + d] = rot_per_node
            end
        end
    end
    
    return elem
end

"""Assign global DOF indices for a CompositeShellTri3 element."""
function populate_globalID!(elem::CompositeShellTri3)
    elem.globalID = vcat([n.globalID[1:6] for n in elem.nodes]...)
    return elem
end

# =============================================================================
# Stress/Moment Recovery
# =============================================================================

"""
    membrane_forces(elem::CompositeShellTri3, u_global::Vector{Float64}) -> Vector{Float64}

Compute membrane forces [Nxx, Nyy, Nxy] in local coordinates (N/m).
"""
function membrane_forces(elem::CompositeShellTri3, u_global::Vector{Float64})
    u_elem = zeros(18)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:6
            global_dof = node.globalID[j]
            u_elem[(i-1)*6 + j] = u_global[global_dof]
        end
    end
    
    u_local = elem.R * u_elem
    
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    Bm = zeros(3, 18)
    Bb = zeros(3, 18)
    _Bmmat!(Bm, gradN_e)
    _Bbmat!(Bb, gradN_e)
    
    ε = Bm * u_local  # Mid-plane strains
    κ = Bb * u_local  # Curvatures
    
    A, B, _ = laminate_stiffnesses(elem.laminate)
    
    # N = A*ε + B*κ
    N = A * ε + B * κ
    
    return N
end

"""
    bending_moments(elem::CompositeShellTri3, u_global::Vector{Float64}) -> Vector{Float64}

Compute bending moments [Mxx, Myy, Mxy] in local coordinates (N·m/m).
"""
function bending_moments(elem::CompositeShellTri3, u_global::Vector{Float64})
    u_elem = zeros(18)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:6
            global_dof = node.globalID[j]
            u_elem[(i-1)*6 + j] = u_global[global_dof]
        end
    end
    
    u_local = elem.R * u_elem
    
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    Bm = zeros(3, 18)
    Bb = zeros(3, 18)
    _Bmmat!(Bm, gradN_e)
    _Bbmat!(Bb, gradN_e)
    
    ε = Bm * u_local
    κ = Bb * u_local
    
    _, B, D = laminate_stiffnesses(elem.laminate)
    
    # M = B*ε + D*κ
    M = B * ε + D * κ
    
    return M
end

"""
    ply_stresses(elem::CompositeShellTri3, u_global::Vector{Float64}, ply_idx::Int, z_position::Symbol=:middle)

Compute stresses in a specific ply at given through-thickness position.

# Arguments
- `ply_idx`: Index of ply (1 = bottom ply)
- `z_position`: `:top`, `:middle`, or `:bottom` of the ply

Returns: [σ1, σ2, τ12] in ply material coordinates (Pa)
"""
function ply_stresses(elem::CompositeShellTri3, u_global::Vector{Float64}, 
                      ply_idx::Int; z_position::Symbol=:middle)
    lam = elem.laminate
    @assert 1 <= ply_idx <= length(lam.plies) "Invalid ply index"
    
    # Get strains and curvatures
    u_elem = zeros(18)
    for (i, node) in enumerate(elem.nodes)
        for j in 1:6
            global_dof = node.globalID[j]
            u_elem[(i-1)*6 + j] = u_global[global_dof]
        end
    end
    u_local = elem.R * u_elem
    
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    Bm = zeros(3, 18)
    Bb = zeros(3, 18)
    _Bmmat!(Bm, gradN_e)
    _Bbmat!(Bb, gradN_e)
    
    ε0 = Bm * u_local  # Mid-plane strains
    κ = Bb * u_local   # Curvatures
    
    # Find z-coordinate of ply
    t_total = thickness(lam)
    z = -t_total / 2 - lam.offset
    for i in 1:(ply_idx-1)
        z += lam.plies[i].thickness
    end
    
    ply = lam.plies[ply_idx]
    z_local = if z_position == :bottom
        z
    elseif z_position == :top
        z + ply.thickness
    else  # :middle
        z + ply.thickness / 2
    end
    
    # Strain at z: ε(z) = ε0 + z*κ
    ε_z = ε0 + z_local * κ
    
    # Transform strain from laminate to ply coordinates
    θ = ply.angle * π / 180
    Tbar = _strain_transform_matrix(θ)
    ε_ply = Tbar * ε_z
    
    # Stress in ply coordinates: σ = Q * ε
    σ_ply = ply.Qbar * ε_ply
    
    return σ_ply
end

# Note: SurfaceLoad works with CompositeShellTri3 automatically since it's a ShellElement
