#=
Stiffness Matrix Module for Asap.jl
===================================

Provides element stiffness matrices using the natural deformation formulation
for improved numerical stability, with both Bernoulli-Euler (shear-rigid) and
Timoshenko (shear-flexible) beam formulations.

Natural Deformation Formulation:
The 12 Cartesian DOFs are transformed to 6 natural deformations:
- dN[1]: Axial elongation
- dN[2]: Symmetric bending in x1-x2 plane
- dN[3]: Anti-symmetric bending in x1-x2 plane  
- dN[4]: Symmetric bending in x1-x3 plane
- dN[5]: Anti-symmetric bending in x1-x3 plane
- dN[6]: Axial torque

This provides better numerical conditioning and clearer physical interpretation.

Eccentric connection support transforms element stiffness to account for
offset connections from node positions.

Mathematical formulations based on:
[1] J.H. Argyris, P.C. Dunne, D.W. Scharpf, "On large displacement-small strain 
    analysis of structures with rotational DOFs", CMAME 14 (1978) 401-451
[2] FinEtoolsFlexStructures.jl by Petr Krysl - MIT License

=#

using LinearAlgebra: mul!, Transpose, I

# =============================================================================
# Natural Deformation Formulation
# =============================================================================

"""
    local_cartesian_to_natural!(aN, L)

Compute transformation matrix from local Cartesian displacements to natural deformations.

Matrix defined in Eq (4.8) of Argyris, Dunne & Scharpf (1978).

# Arguments
- `aN`: Output 6×12 transformation matrix (modified in place)
- `L`: Element length [m]

# Natural deformations
1. Axial elongation
2. Symmetric bending (x1-x2)
3. Anti-symmetric bending (x1-x2)
4. Symmetric bending (x1-x3)
5. Anti-symmetric bending (x1-x3)
6. Torsional twist
"""
function local_cartesian_to_natural!(aN::Matrix{Float64}, L)
    fill!(aN, 0.0)
    
    aN[1, 1] = -1;  aN[1, 7] = +1                           # Axial
    aN[2, 6] = +1;  aN[2, 12] = -1                          # Symmetric bending (x1-x2)
    aN[3, 2] = 2/L; aN[3, 6] = +1; aN[3, 8] = -2/L; aN[3, 12] = +1  # Anti-sym bending (x1-x2)
    aN[4, 5] = -1;  aN[4, 11] = +1                          # Symmetric bending (x1-x3)
    aN[5, 3] = 2/L; aN[5, 5] = -1; aN[5, 9] = -2/L; aN[5, 11] = -1  # Anti-sym bending (x1-x3)
    aN[6, 4] = -1;  aN[6, 10] = +1                          # Torsion
    
    return aN
end

"""
    natural_stiffness_bernoulli!(DN, E, G, A, I2, I3, J, L)

Natural stiffness matrix for Bernoulli-Euler (shear-rigid) beam.

# Arguments
- `DN`: Output 6×6 natural stiffness matrix
- `E`: Young's modulus [Pa]
- `G`: Shear modulus [Pa]  
- `A`: Cross-sectional area [m²]
- `I2`: Weak axis moment of inertia [m⁴]
- `I3`: Strong axis moment of inertia [m⁴]
- `J`: Torsional constant [m⁴]
- `L`: Element length [m]
"""
function natural_stiffness_bernoulli!(DN::Matrix{Float64}, E, G, A, I2, I3, J, L)
    fill!(DN, 0.0)
    DN[1, 1] = E * A / L        # Axial
    DN[2, 2] = E * I3 / L       # Symmetric bending (strong)
    DN[3, 3] = 3 * E * I3 / L   # Anti-symmetric bending (strong)
    DN[4, 4] = E * I2 / L       # Symmetric bending (weak)
    DN[5, 5] = 3 * E * I2 / L   # Anti-symmetric bending (weak)
    DN[6, 6] = G * J / L        # Torsion
    return DN
end

"""
    natural_stiffness_timoshenko!(DN, E, G, A, I2, I3, J, Ay, Az, L)

Natural stiffness matrix for Timoshenko (shear-flexible) beam.

The shear correction modifies the anti-symmetric bending stiffnesses.

# Arguments
- `DN`: Output 6×6 natural stiffness matrix
- `E`: Young's modulus [Pa]
- `G`: Shear modulus [Pa]
- `A`: Cross-sectional area [m²]  
- `I2`: Weak axis moment of inertia [m⁴]
- `I3`: Strong axis moment of inertia [m⁴]
- `J`: Torsional constant [m⁴]
- `Ay`: Shear area in y-direction [m²]
- `Az`: Shear area in z-direction [m²]
- `L`: Element length [m]
"""
function natural_stiffness_timoshenko!(DN::Matrix{Float64}, E, G, A, I2, I3, J, Ay, Az, L)
    fill!(DN, 0.0)
    
    DN[1, 1] = E * A / L        # Axial
    DN[2, 2] = E * I3 / L       # Symmetric bending (strong)
    DN[4, 4] = E * I2 / L       # Symmetric bending (weak)
    DN[6, 6] = G * J / L        # Torsion
    
    # Shear-corrected anti-symmetric bending
    # Phi = 12*E*I / (G*A_shear*L²) is the shear correction factor
    Phi3 = 12 * E * I3 / (G * Ay * L^2)
    Phi2 = 12 * E * I2 / (G * Az * L^2)
    
    DN[3, 3] = 3 * E * I3 / L / (1 + Phi3)  # Anti-sym (strong)
    DN[5, 5] = 3 * E * I2 / L / (1 + Phi2)  # Anti-sym (weak)
    
    return DN
end

"""
    natural_stiffness!(DN, E, G, A, I2, I3, J, Ay, Az, L)

Compute natural stiffness matrix, automatically selecting Bernoulli-Euler
or Timoshenko formulation based on shear areas.
"""
function natural_stiffness!(DN::Matrix{Float64}, E, G, A, I2, I3, J, Ay, Az, L)
    if !isfinite(Ay) || !isfinite(Az)
        return natural_stiffness_bernoulli!(DN, E, G, A, I2, I3, J, L)
    else
        return natural_stiffness_timoshenko!(DN, E, G, A, I2, I3, J, Ay, Az, L)
    end
end

"""
    local_stiffness_natural!(K, E, G, A, I2, I3, J, Ay, Az, L, aN, DN)

Compute local stiffness matrix using natural deformation formulation.

K = aN' * DN * aN

# Arguments
- `K`: Output 12×12 stiffness matrix
- Other parameters: Material and geometric properties
- `aN`: Pre-allocated 6×12 transformation matrix
- `DN`: Pre-allocated 6×6 natural stiffness matrix
"""
function local_stiffness_natural!(K::Matrix{Float64}, E, G, A, I2, I3, J, Ay, Az, L,
                                  aN::Matrix{Float64}, DN::Matrix{Float64})
    local_cartesian_to_natural!(aN, L)
    natural_stiffness!(DN, E, G, A, I2, I3, J, Ay, Az, L)
    K .= aN' * DN * aN
    return K
end

# =============================================================================
# Eccentricity Transformation
# =============================================================================

"""
    eccentricity_transformation!(Te, e_f1_1, e_f1_2, e_f2, e_f3)

Compute transformation matrix for eccentric beam connections.

When a beam is offset from the node positions, the stiffness must be
transformed to account for the rigid link between the node and beam end.

# Arguments
- `Te`: Output 12×12 transformation matrix
- `e_f1_1`: Axial eccentricity at start node [m]
- `e_f1_2`: Axial eccentricity at end node [m]
- `e_f2`: Transverse eccentricity in f2 (y) direction [m]
- `e_f3`: Transverse eccentricity in f3 (z) direction [m]

The transformation relates master (node) DOFs to slave (beam end) DOFs:
U_slave = Te * U_master
"""
function eccentricity_transformation!(Te::Matrix{Float64}, e_f1_1, e_f1_2, e_f2, e_f3)
    Te .= 0.0
    
    # Identity blocks for direct DOF connections
    Te[1:3, 1:3] .= I(3)
    Te[4:6, 4:6] .= I(3)
    Te[7:9, 7:9] .= I(3)
    Te[10:12, 10:12] .= I(3)
    
    # Coupling terms from rotation to translation (rigid body kinematics)
    # At start node (rows 1-3 from cols 4-6)
    Te[1, 4:6] .= (0.0, e_f3, -e_f2)
    Te[2, 4:6] .= (-e_f3, 0.0, e_f1_1)
    Te[3, 4:6] .= (e_f2, -e_f1_1, 0.0)
    
    # At end node (rows 7-9 from cols 10-12)
    Te[7, 10:12] .= (0.0, e_f3, -e_f2)
    Te[8, 10:12] .= (-e_f3, 0.0, e_f1_2)
    Te[9, 10:12] .= (e_f2, -e_f1_2, 0.0)
    
    return Te
end

# =============================================================================
# Direct Stiffness Matrix Functions (Original ASAP API - kept for compatibility)
# =============================================================================

"""
    k_fixedfixed(E, A, L, G, Ix, Iy, J; Ay=Inf, Az=Inf)

12×12 stiffness matrix for beam with full coupling at both ends.
Uses Timoshenko formulation when Ay, Az are finite.
"""
function k_fixedfixed(E, A, L, G, Ix, Iy, J; Ay=Inf, Az=Inf)
    # Use natural formulation
    aN = zeros(6, 12)
    DN = zeros(6, 6)
    K = zeros(12, 12)
    local_stiffness_natural!(K, E, G, A, Iy, Ix, J, Ay, Az, L, aN, DN)
    return K
end

"""
    k_freefixed(E, A, L, Ix, Iy)

12×12 stiffness matrix with rotational DOFs released at start node.
"""
function k_freefixed(E, A, L, Ix, Iy)
    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 0 0 -3Ix 0 0 0 3L*Ix;
        0 0 3Iy 0 0 0 0 0 -3Iy 0 -3L*Iy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 0 0 3Ix 0 0 0 -3L*Ix;
        0 0 -3Iy 0 0 0 0 0 3Iy 0 3L*Iy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 0 0 0 0 3L*Iy 0 3L^2*Iy 0;
        0 3L*Ix 0 0 0 0 0 -3L*Ix 0 0 0 3L^2*Ix    
    ]
    return k
end

"""
    k_fixedfree(E, A, L, Ix, Iy)

12×12 stiffness matrix with rotational DOFs released at end node.
"""
function k_fixedfree(E, A, L, Ix, Iy)
    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 3L*Ix 0 -3Ix 0 0 0 0;
        0 0 3Iy 0 -3L*Iy 0 0 0 -3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 3L^2*Iy 0 0 0 3L*Iy 0 0 0;
        0 3L*Ix 0 0 0 3L^2*Ix 0 -3L*Ix 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 -3L*Ix 0 3Ix 0 0 0 0;
        0 0 -3Iy 0 3L*Iy 0 0 0 3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
    ]
    return k
end

"""
    k_freefree(E, A, L)

12×12 stiffness matrix for truss-like behavior (axial only).
"""
function k_freefree(E, A, L)
    k = zeros(12, 12)
    k[1,1] = 1
    k[1,7] = -1
    k[7,1] = -1
    k[7,7] = 1
    return E * A / L .* k
end

"""
    k_joist(E, A, L, G, J)

12×12 stiffness matrix with only axial and torsional coupling.
"""
function k_joist(E, A, L, G, J)
    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
    ]
    return k
end

# =============================================================================
# Element Local Stiffness Functions
# =============================================================================

"""
    local_K(element::Element)

Return the element stiffness matrix in local coordinate system.
"""
local_K(element::Element{FixedFixed}) = k_fixedfixed(
    to_pascals(element.section.E),
    to_meters_squared(element.section.A),
    to_meters(element.length),
    to_pascals(element.section.G),
    to_meters_fourth(element.section.Ix),
    to_meters_fourth(element.section.Iy),
    to_meters_fourth(element.section.J);
    Ay = element.section.Ay,
    Az = element.section.Az
)

local_K(element::Element{FixedFree}) = k_fixedfree(
    to_pascals(element.section.E),
    to_meters_squared(element.section.A),
    to_meters(element.length),
    to_meters_fourth(element.section.Ix),
    to_meters_fourth(element.section.Iy)
)

local_K(element::Element{FreeFixed}) = k_freefixed(
    to_pascals(element.section.E),
    to_meters_squared(element.section.A),
    to_meters(element.length),
    to_meters_fourth(element.section.Ix),
    to_meters_fourth(element.section.Iy)
)

local_K(element::Element{FreeFree}) = k_freefree(
    to_pascals(element.section.E),
    to_meters_squared(element.section.A),
    to_meters(element.length)
)

local_K(element::Element{Joist}) = k_joist(
    to_pascals(element.section.E),
    to_meters_squared(element.section.A),
    to_meters(element.length),
    to_pascals(element.section.G),
    to_meters_fourth(element.section.J)
)

local_K(element::TrussElement) = to_pascals(element.section.E) * to_meters_squared(element.section.A) / to_meters(element.length) .* [1 -1; -1 1]

# =============================================================================
# Global Stiffness Matrix Functions (with eccentricity support)
# =============================================================================

"""
    global_K!(element::Element)

Populate the element stiffness matrix `element.K` in global coordinate system.
Applies eccentricity transformation if element has eccentric connections.
"""
function global_K!(element::Element)
    K_local = local_K(element)
    
    # Apply eccentricity transformation if needed
    if has_eccentricity(element)
        Te = zeros(12, 12)
        eccentricity_transformation!(Te, element.eccentricity...)
        K_temp = K_local * Te
        K_local = Te' * K_temp
    end
    
    # Transform to global coordinates
    element.K = element.R' * K_local * element.R
end

function global_K!(element::TrussElement)
    element.K = element.R' * local_K(element) * element.R
end

"""
    global_K(element::Element)

Return the element stiffness matrix in global coordinate system.
"""
function global_K(element::Element)
    K_local = local_K(element)
    
    if has_eccentricity(element)
        Te = zeros(12, 12)
        eccentricity_transformation!(Te, element.eccentricity...)
        K_local = Te' * K_local * Te
    end
    
    return element.R' * K_local * element.R
end

function global_K(element::TrussElement)
    return element.R' * local_K(element) * element.R
end
