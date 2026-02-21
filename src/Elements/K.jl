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
# Thread-Local Workspace Buffers
# =============================================================================

struct _KWorkspace
    aN::Matrix{Float64}       # 6×12
    DN::Matrix{Float64}       # 6×6
    K_local::Matrix{Float64}  # 12×12
    temp6x12::Matrix{Float64} # 6×12  (for mul! in natural stiffness)
    temp12x12::Matrix{Float64}# 12×12 (for mul! in global_K!)
    Te::Matrix{Float64}       # 12×12 (for eccentricity transformation)
end

_KWorkspace() = _KWorkspace(
    zeros(6, 12), zeros(6, 6), zeros(12, 12),
    zeros(6, 12), zeros(12, 12), zeros(12, 12)
)

const _K_WS_POOL = [_KWorkspace() for _ in 1:max(1, Threads.nthreads())]

@inline function _get_k_ws()
    tid = Threads.threadid()
    @inbounds return _K_WS_POOL[tid]
end

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
    local_stiffness_natural!(K, E, G, A, I2, I3, J, Ay, Az, L, aN, DN, temp6x12)

Compute local stiffness matrix using natural deformation formulation.
Uses mul! to avoid intermediate allocations:
  K = aN' * DN * aN  →  temp = DN * aN;  K = aN' * temp

# Arguments
- `K`: Output 12×12 stiffness matrix
- Other parameters: Material and geometric properties
- `aN`: Pre-allocated 6×12 transformation matrix
- `DN`: Pre-allocated 6×6 natural stiffness matrix
- `temp6x12`: Pre-allocated 6×12 workspace
"""
function local_stiffness_natural!(K::Matrix{Float64}, E, G, A, I2, I3, J, Ay, Az, L,
                                  aN::Matrix{Float64}, DN::Matrix{Float64},
                                  temp6x12::Matrix{Float64})
    local_cartesian_to_natural!(aN, L)
    natural_stiffness!(DN, E, G, A, I2, I3, J, Ay, Az, L)
    mul!(temp6x12, DN, aN)          # 6×12 = 6×6 × 6×12
    mul!(K, transpose(aN), temp6x12) # 12×12 = 12×6 × 6×12
    return K
end

# Overload without pre-allocated temp buffer (allocates temp6x12 internally)
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
    @inbounds begin
        Te[1,1] = 1.0; Te[2,2] = 1.0; Te[3,3] = 1.0
        Te[4,4] = 1.0; Te[5,5] = 1.0; Te[6,6] = 1.0
        Te[7,7] = 1.0; Te[8,8] = 1.0; Te[9,9] = 1.0
        Te[10,10] = 1.0; Te[11,11] = 1.0; Te[12,12] = 1.0
        
        # Coupling at start node (rows 1-3 from cols 4-6)
        Te[1, 4] = 0.0;    Te[1, 5] = e_f3;    Te[1, 6] = -e_f2
        Te[2, 4] = -e_f3;  Te[2, 5] = 0.0;     Te[2, 6] = e_f1_1
        Te[3, 4] = e_f2;   Te[3, 5] = -e_f1_1; Te[3, 6] = 0.0
        
        # Coupling at end node (rows 7-9 from cols 10-12)
        Te[7, 10] = 0.0;    Te[7, 11] = e_f3;    Te[7, 12] = -e_f2
        Te[8, 10] = -e_f3;  Te[8, 11] = 0.0;     Te[8, 12] = e_f1_2
        Te[9, 10] = e_f2;   Te[9, 11] = -e_f1_2; Te[9, 12] = 0.0
    end
    
    return Te
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
    Ay = to_meters_squared(element.section.Ay),
    Az = to_meters_squared(element.section.Az)
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
    global_K!(element::Element{FixedFixed})

Fully allocation-free stiffness computation for the most common release type.
Computes local K via natural deformation, applies eccentricity if needed,
then transforms to global coordinates using mul! with thread-local workspace.
"""
function global_K!(element::Element{FixedFixed})
    ws = _get_k_ws()
    
    E  = to_pascals(element.section.E)
    A  = to_meters_squared(element.section.A)
    L  = to_meters(element.length)
    G  = to_pascals(element.section.G)
    Ix = to_meters_fourth(element.section.Ix)
    Iy = to_meters_fourth(element.section.Iy)
    J  = to_meters_fourth(element.section.J)
    
    local_stiffness_natural!(ws.K_local, E, G, A, Iy, Ix, J,
                             to_meters_squared(element.section.Ay), to_meters_squared(element.section.Az), L,
                             ws.aN, ws.DN, ws.temp6x12)
    
    if has_eccentricity(element)
        eccentricity_transformation!(ws.Te, element.eccentricity...)
        mul!(ws.temp12x12, ws.K_local, ws.Te)
        mul!(ws.K_local, transpose(ws.Te), ws.temp12x12)
    end
    
    mul!(ws.temp12x12, ws.K_local, element.R)
    mul!(element.K, transpose(element.R), ws.temp12x12)
end

"""
    global_K!(element::Element)

Generic global stiffness for non-FixedFixed release types.
Uses mul! for the R' * K_local * R transformation to reduce allocations.
"""
function global_K!(element::Element)
    ws = _get_k_ws()
    K_local = local_K(element)
    
    if has_eccentricity(element)
        eccentricity_transformation!(ws.Te, element.eccentricity...)
        mul!(ws.temp12x12, K_local, ws.Te)
        copyto!(ws.K_local, ws.temp12x12)
        # K_local is now ws.temp12x12, but we need Te' * that
        # Use K_local as scratch
        mul!(ws.temp12x12, transpose(ws.Te), ws.K_local)
        K_local = ws.temp12x12
    end
    
    # Avoid using ws.temp12x12 if K_local IS ws.temp12x12
    # Safe path: copy into ws.K_local first
    if K_local === ws.temp12x12
        copyto!(ws.K_local, K_local)
        mul!(ws.temp12x12, ws.K_local, element.R)
    else
        mul!(ws.temp12x12, K_local, element.R)
    end
    mul!(element.K, transpose(element.R), ws.temp12x12)
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
