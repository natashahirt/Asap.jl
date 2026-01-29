# =============================================================================
# Shell Elements for Membrane/Diaphragm Analysis
# =============================================================================

abstract type ShellElement <: AbstractElement end

# =============================================================================
# ShellTri3: Constant Strain Triangle (CST) Membrane Element
# =============================================================================

"""
    ShellTri3(nodes, thickness, E, ν, id=:shell)

3-node triangular membrane element (constant strain triangle).
In-plane behavior only (2 DOF per node in element plane → 6 DOF total).
Uses global DOFs 1,2,3 (translations) at each node.

# Arguments
- `nodes`: Tuple of 3 Node objects defining the triangle
- `thickness`: Element thickness (Unitful length)
- `E`: Young's modulus (Unitful pressure)
- `ν`: Poisson's ratio
"""
mutable struct ShellTri3 <: ShellElement
    nodes::NTuple{3, Node}
    thickness::QuantityDistance
    E::QuantityPressure
    ν::Float64
    elementID::Int64
    globalID::Vector{Int64}
    K::Matrix{Float64}      # 9×9 in GCS (3 translations per node)
    R::Matrix{Float64}      # 6×9 transformation (local 2D → global 3D)
    LCS::NTuple{3, Vector{Float64}}  # Local coordinate system (x, y, normal)
    area::Float64           # Element area (m²)
    id::Symbol

    function ShellTri3(nodes::NTuple{3, Node}, thickness, E, ν, id=:shell)
        @assert 0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        
        elem = new(
            nodes,
            uconvert(u"m", thickness),
            uconvert(u"Pa", E),
            ν,
            0,
            Vector{Int64}(undef, 9),
            zeros(9, 9),
            zeros(6, 9),
            (zeros(3), zeros(3), zeros(3)),
            0.0,
            id
        )
        return elem
    end
end

# Convenience constructor from vector
ShellTri3(nodes::Vector{Node}, thickness, E, ν, id=:shell) = 
    ShellTri3((nodes[1], nodes[2], nodes[3]), thickness, E, ν, id)

# =============================================================================
# Geometry & Local Coordinate System
# =============================================================================

"""Compute local coordinate system for a shell element from its nodes."""
function lcs!(elem::ShellTri3)
    p1 = [to_meters(c) for c in elem.nodes[1].position]
    p2 = [to_meters(c) for c in elem.nodes[2].position]
    p3 = [to_meters(c) for c in elem.nodes[3].position]
    
    # Local x: along edge 1→2
    v12 = p2 - p1
    x_local = normalize(v12)
    
    # Edge 1→3
    v13 = p3 - p1
    
    # Normal (local z): cross product
    n = cross(v12, v13)
    area = norm(n) / 2
    z_local = normalize(n)
    
    # Local y: perpendicular to x in the plane
    y_local = cross(z_local, x_local)
    
    elem.LCS = (x_local, y_local, z_local)
    elem.area = area
end

"""Compute element area."""
function area(elem::ShellTri3)
    if elem.area == 0.0
        lcs!(elem)
    end
    return elem.area
end

# =============================================================================
# Transformation Matrix
# =============================================================================

"""
Build transformation matrix R for ShellTri3.
Maps 6 local DOFs (2 per node in-plane) to 9 global DOFs (3 per node in 3D).

Local DOFs: [u1_x, u1_y, u2_x, u2_y, u3_x, u3_y] (in element plane)
Global DOFs: [u1_X, u1_Y, u1_Z, u2_X, u2_Y, u2_Z, u3_X, u3_Y, u3_Z]
"""
function R!(elem::ShellTri3)
    x_loc, y_loc, _ = elem.LCS
    
    # 2×3 block: transforms global 3D displacement to local 2D in-plane
    # [u_local_x]   [x_loc']   [u_global]
    # [u_local_y] = [y_loc'] × [v_global]
    #                          [w_global]
    T_block = [x_loc'; y_loc']  # 2×3
    
    # Full transformation: 6×9 (3 nodes × 2 local DOF, 3 nodes × 3 global DOF)
    elem.R .= 0.0
    for i in 1:3
        rows = (2*(i-1)+1):(2*i)
        cols = (3*(i-1)+1):(3*i)
        elem.R[rows, cols] = T_block
    end
end

# =============================================================================
# Stiffness Matrix
# =============================================================================

"""
Compute local stiffness matrix for ShellTri3 (6×6 in local coords).
Constant Strain Triangle formulation.
"""
function local_K(elem::ShellTri3)
    # Node coordinates in local system
    p1 = [to_meters(c) for c in elem.nodes[1].position]
    p2 = [to_meters(c) for c in elem.nodes[2].position]
    p3 = [to_meters(c) for c in elem.nodes[3].position]
    
    x_loc, y_loc, _ = elem.LCS
    
    # Project to local 2D coordinates
    x = [dot(p - p1, x_loc) for p in [p1, p2, p3]]
    y = [dot(p - p1, y_loc) for p in [p1, p2, p3]]
    
    # Area (should match elem.area)
    A = elem.area
    
    # B matrix coefficients (strain-displacement)
    # β_i = y_j - y_k, γ_i = x_k - x_j (cyclic)
    β = [y[2] - y[3], y[3] - y[1], y[1] - y[2]]
    γ = [x[3] - x[2], x[1] - x[3], x[2] - x[1]]
    
    # B matrix: ε = B × u_local
    # [εxx]       1   [β1  0  β2  0  β3  0 ]   [u1]
    # [εyy] = ―――――  [ 0 γ1  0  γ2  0  γ3]   [v1]
    # [γxy]     2A   [γ1 β1 γ2 β2 γ3 β3]   [...]
    B = (1 / (2A)) * [
        β[1]  0     β[2]  0     β[3]  0    ;
        0     γ[1]  0     γ[2]  0     γ[3] ;
        γ[1]  β[1]  γ[2]  β[2]  γ[3]  β[3]
    ]
    
    # Constitutive matrix (plane stress)
    t = to_meters(elem.thickness)
    E = to_pascals(elem.E)
    ν = elem.ν
    
    D = (E / (1 - ν^2)) * [
        1  ν  0       ;
        ν  1  0       ;
        0  0  (1-ν)/2
    ]
    
    # Stiffness: K = t × A × B' × D × B
    K_local = t * A * (B' * D * B)
    
    return K_local
end

"""Populate global stiffness matrix K for ShellTri3."""
function global_K!(elem::ShellTri3)
    K_local = local_K(elem)  # 6×6
    # K_global = R' × K_local × R  (9×6 × 6×6 × 6×9 = 9×9)
    elem.K = elem.R' * K_local * elem.R
end

# =============================================================================
# Processing Interface (matches frame elements)
# =============================================================================

"""Process a shell element: compute LCS, R, and K."""
function process!(elem::ShellTri3)
    lcs!(elem)
    R!(elem)
    global_K!(elem)
end

# =============================================================================
# ShellQuad4: 4-Node Isoparametric Membrane Element
# =============================================================================

"""
    ShellQuad4(nodes, thickness, E, ν, id=:shell)

4-node quadrilateral membrane element (bilinear isoparametric).
In-plane behavior only (2 DOF per node in element plane → 8 DOF total).
Uses global DOFs 1,2,3 (translations) at each node.

# Arguments
- `nodes`: Tuple of 4 Node objects defining the quad (counter-clockwise)
- `thickness`: Element thickness (Unitful length)
- `E`: Young's modulus (Unitful pressure)
- `ν`: Poisson's ratio
"""
mutable struct ShellQuad4 <: ShellElement
    nodes::NTuple{4, Node}
    thickness::QuantityDistance
    E::QuantityPressure
    ν::Float64
    elementID::Int64
    globalID::Vector{Int64}
    K::Matrix{Float64}      # 12×12 in GCS (3 translations per node)
    R::Matrix{Float64}      # 8×12 transformation (local 2D → global 3D)
    LCS::NTuple{3, Vector{Float64}}  # Local coordinate system (x, y, normal)
    area::Float64           # Element area (m²)
    id::Symbol

    function ShellQuad4(nodes::NTuple{4, Node}, thickness, E, ν, id=:shell)
        @assert 0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        
        elem = new(
            nodes,
            uconvert(u"m", thickness),
            uconvert(u"Pa", E),
            ν,
            0,
            Vector{Int64}(undef, 12),
            zeros(12, 12),
            zeros(8, 12),
            (zeros(3), zeros(3), zeros(3)),
            0.0,
            id
        )
        return elem
    end
end

# Convenience constructor
ShellQuad4(nodes::Vector{Node}, thickness, E, ν, id=:shell) = 
    ShellQuad4((nodes[1], nodes[2], nodes[3], nodes[4]), thickness, E, ν, id)

"""Compute local coordinate system for ShellQuad4 from centroid and edge vectors."""
function lcs!(elem::ShellQuad4)
    p = [[to_meters(c) for c in n.position] for n in elem.nodes]
    
    # Approximate plane normal from cross products of diagonals
    d13 = p[3] - p[1]
    d24 = p[4] - p[2]
    n = cross(d13, d24)
    z_local = normalize(n)
    
    # Local x: average of bottom and top edges
    v12 = p[2] - p[1]
    v43 = p[3] - p[4]
    x_local = normalize(v12 + v43)
    
    # Local y: perpendicular to x in the plane
    y_local = cross(z_local, x_local)
    
    elem.LCS = (x_local, y_local, z_local)
    
    # Compute area using shoelace formula in local coords
    x_loc = [dot(pi - p[1], x_local) for pi in p]
    y_loc = [dot(pi - p[1], y_local) for pi in p]
    elem.area = 0.5 * abs(
        (x_loc[1] - x_loc[3]) * (y_loc[2] - y_loc[4]) - 
        (x_loc[2] - x_loc[4]) * (y_loc[1] - y_loc[3])
    )
end

"""Build transformation matrix R for ShellQuad4 (8×12)."""
function R!(elem::ShellQuad4)
    x_loc, y_loc, _ = elem.LCS
    T_block = [x_loc'; y_loc']  # 2×3
    
    elem.R .= 0.0
    for i in 1:4
        rows = (2*(i-1)+1):(2*i)
        cols = (3*(i-1)+1):(3*i)
        elem.R[rows, cols] = T_block
    end
end

"""
Bilinear shape functions and derivatives at natural coordinates (ξ, η).
Returns (N, dNdξ, dNdη) where N is 4-element, dNdξ/dNdη are 4-element.
"""
function shape_functions_quad4(ξ, η)
    # Shape functions at corners: N_i = (1/4)(1 + ξ_i*ξ)(1 + η_i*η)
    # Corner natural coords: (-1,-1), (1,-1), (1,1), (-1,1)
    ξ_c = [-1, 1, 1, -1]
    η_c = [-1, -1, 1, 1]
    
    N = [(1 + ξ_c[i]*ξ) * (1 + η_c[i]*η) / 4 for i in 1:4]
    dNdξ = [ξ_c[i] * (1 + η_c[i]*η) / 4 for i in 1:4]
    dNdη = [(1 + ξ_c[i]*ξ) * η_c[i] / 4 for i in 1:4]
    
    return N, dNdξ, dNdη
end

"""Compute local stiffness matrix for ShellQuad4 (8×8) using 2×2 Gauss quadrature."""
function local_K(elem::ShellQuad4)
    # Node coordinates in local 2D system
    p = [[to_meters(c) for c in n.position] for n in elem.nodes]
    x_loc, y_loc, _ = elem.LCS
    
    x = [dot(pi - p[1], x_loc) for pi in p]
    y = [dot(pi - p[1], y_loc) for pi in p]
    
    # Material
    t = to_meters(elem.thickness)
    E = to_pascals(elem.E)
    ν = elem.ν
    
    # Plane stress constitutive matrix
    D = (E / (1 - ν^2)) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    
    # 2×2 Gauss quadrature
    gp = 1 / sqrt(3)
    gauss_pts = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
    
    K = zeros(8, 8)
    
    for (ξ, η) in gauss_pts
        N, dNdξ, dNdη = shape_functions_quad4(ξ, η)
        
        # Jacobian matrix
        # J = [∂x/∂ξ  ∂y/∂ξ; ∂x/∂η  ∂y/∂η]
        dxdξ = sum(dNdξ .* x)
        dydξ = sum(dNdξ .* y)
        dxdη = sum(dNdη .* x)
        dydη = sum(dNdη .* y)
        
        J = [dxdξ dydξ; dxdη dydη]
        detJ = det(J)
        Jinv = inv(J)
        
        # Shape function derivatives in physical coords
        # [dN/dx; dN/dy] = J^(-1) * [dN/dξ; dN/dη]
        dNdx = Jinv[1,1] .* dNdξ + Jinv[1,2] .* dNdη
        dNdy = Jinv[2,1] .* dNdξ + Jinv[2,2] .* dNdη
        
        # B matrix (3×8): strain-displacement
        # [εxx]   [dN1/dx   0    dN2/dx   0    ...]
        # [εyy] = [  0   dN1/dy    0   dN2/dy  ...] × [u1, v1, u2, v2, ...]
        # [γxy]   [dN1/dy dN1/dx dN2/dy dN2/dx ...]
        B = zeros(3, 8)
        for i in 1:4
            B[1, 2i-1] = dNdx[i]
            B[2, 2i]   = dNdy[i]
            B[3, 2i-1] = dNdy[i]
            B[3, 2i]   = dNdx[i]
        end
        
        # Integrate: K += B' * D * B * t * detJ * weight
        K += t * detJ * (B' * D * B)  # weight = 1 for 2×2 Gauss
    end
    
    return K
end

"""Populate global stiffness matrix K for ShellQuad4."""
function global_K!(elem::ShellQuad4)
    K_local = local_K(elem)  # 8×8
    elem.K = elem.R' * K_local * elem.R  # 12×12
end

"""Process ShellQuad4: compute LCS, R, and K."""
function process!(elem::ShellQuad4)
    lcs!(elem)
    R!(elem)
    global_K!(elem)
end

# =============================================================================
# Stress Recovery
# =============================================================================

"""
Compute membrane stresses for ShellQuad4 at element centroid.
Returns [σxx, σyy, τxy] in Pa.
"""
function stress(elem::ShellQuad4, u_global::Vector{Float64})
    idx = elem.globalID
    u_elem = u_global[idx]  # 12 global DOFs
    u_local = elem.R * u_elem  # 8 local DOFs
    
    # Get local coordinates
    p = [[to_meters(c) for c in n.position] for n in elem.nodes]
    x_loc, y_loc, _ = elem.LCS
    x = [dot(pi - p[1], x_loc) for pi in p]
    y = [dot(pi - p[1], y_loc) for pi in p]
    
    # Evaluate at centroid (ξ=0, η=0)
    _, dNdξ, dNdη = shape_functions_quad4(0.0, 0.0)
    
    dxdξ = sum(dNdξ .* x)
    dydξ = sum(dNdξ .* y)
    dxdη = sum(dNdη .* x)
    dydη = sum(dNdη .* y)
    
    J = [dxdξ dydξ; dxdη dydη]
    Jinv = inv(J)
    
    dNdx = Jinv[1,1] .* dNdξ + Jinv[1,2] .* dNdη
    dNdy = Jinv[2,1] .* dNdξ + Jinv[2,2] .* dNdη
    
    B = zeros(3, 8)
    for i in 1:4
        B[1, 2i-1] = dNdx[i]
        B[2, 2i]   = dNdy[i]
        B[3, 2i-1] = dNdy[i]
        B[3, 2i]   = dNdx[i]
    end
    
    ε = B * u_local
    
    E = to_pascals(elem.E)
    ν = elem.ν
    D = (E / (1 - ν^2)) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    
    return D * ε
end

"""
Compute membrane stresses for ShellTri3 from nodal displacements.
Returns [σxx, σyy, τxy] in Pa.
"""
function stress(elem::ShellTri3, u_global::Vector{Float64})
    # Extract element DOFs from global displacement
    idx = elem.globalID
    u_elem = u_global[idx]  # 9 global DOFs
    
    # Transform to local: u_local = R × u_global
    u_local = elem.R * u_elem  # 6 local DOFs
    
    # Compute strains via B matrix
    p1 = [to_meters(c) for c in elem.nodes[1].position]
    p2 = [to_meters(c) for c in elem.nodes[2].position]
    p3 = [to_meters(c) for c in elem.nodes[3].position]
    
    x_loc, y_loc, _ = elem.LCS
    x = [dot(p - p1, x_loc) for p in [p1, p2, p3]]
    y = [dot(p - p1, y_loc) for p in [p1, p2, p3]]
    A = elem.area
    
    β = [y[2] - y[3], y[3] - y[1], y[1] - y[2]]
    γ = [x[3] - x[2], x[1] - x[3], x[2] - x[1]]
    
    B = (1 / (2A)) * [
        β[1]  0     β[2]  0     β[3]  0    ;
        0     γ[1]  0     γ[2]  0     γ[3] ;
        γ[1]  β[1]  γ[2]  β[2]  γ[3]  β[3]
    ]
    
    ε = B * u_local  # [εxx, εyy, γxy]
    
    # Stress via constitutive matrix
    E = to_pascals(elem.E)
    ν = elem.ν
    D = (E / (1 - ν^2)) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    
    σ = D * ε  # [σxx, σyy, τxy] in Pa
    return σ
end

