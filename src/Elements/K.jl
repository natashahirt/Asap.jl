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

# Per-thread pool, filled in `Asap.__init__` (see threading.jl).
const _K_WS_POOL = _KWorkspace[]

@inline _get_k_ws() = _thread_buffer(_K_WS_POOL)

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
# Released-Element Local Stiffness (member end releases)
# =============================================================================
#
# The released-element stiffness matrices below are obtained by *static
# condensation* (member end-release modification) of the fully fixed-fixed
# 12x12 local stiffness matrix. This is a standard, closed-form matrix-analysis
# operation — NOT an approximation and NOT invented mechanics.
#
# Derivation (static condensation of released DOFs):
#   Partition the fixed-fixed stiffness K and the local end forces/displacements
#   into retained (r) and released (c) sets:
#       [f_r]   [K_rr  K_rc] [u_r]
#       [f_c] = [K_cr  K_cc] [u_c]
#   A member end release (internal hinge / moment release) enforces zero force
#   conjugate to each released DOF, i.e. f_c = 0. Solving the second row for u_c
#   and back-substituting gives the reduced stiffness on the retained DOFs:
#       u_c = -K_cc^{-1} K_cr u_r
#       f_r = (K_rr - K_rc K_cc^{-1} K_cr) u_r
#   The released DOFs then carry zero rows and zero columns.
#
# References (member end releases / static condensation):
#   - W. McGuire, R.H. Gallagher, R.D. Ziemian, "Matrix Structural Analysis,"
#     2nd ed., Sec. 7.4 (static condensation) and Sec. 13.x (member releases).
#   - A. Kassimali, "Structural Analysis," 5th ed., Table 8.1 / Table 8.2
#     (member stiffness coefficients with end moment releases).
#   - Weaver & Gere, "Matrix Analysis of Framed Structures," 3rd ed., Sec. 6.5.
#
# Local DOF ordering (matches `local_cartesian_to_natural!` and the assembled
# global DOFs, node i then node j):
#   1: u_x,i  2: u_y,i  3: u_z,i  4: theta_x,i  5: theta_y,i  6: theta_z,i
#   7: u_x,j  8: u_y,j  9: u_z,j 10: theta_x,j 11: theta_y,j 12: theta_z,j
#
# The four physical actions (axial, torsion, bending in the x-y plane, bending
# in the x-z plane) are uncoupled in the local frame, so releasing a rotational
# DOF only affects its own action. Per `natural_stiffness_bernoulli!` and
# `global_K!(::Element{FixedFixed})` the plane assignment is I3 = Ix (bending
# with transverse translation u_y / rotation theta_z) and I2 = Iy (bending with
# transverse translation u_z / rotation theta_y).

"""
    _k_fixedfixed_full(E, A, L, G, Ix, Iy, J, Ay, Az)

Assemble the full 12x12 fixed-fixed local stiffness matrix using the natural
deformation formulation, matching the DOF ordering and I2=Iy / I3=Ix plane
assignment used by `global_K!(::Element{FixedFixed})`.
"""
function _k_fixedfixed_full(E, A, L, G, Ix, Iy, J, Ay, Az)
    K = zeros(12, 12)
    aN = zeros(6, 12)
    DN = zeros(6, 6)
    # I2 = Iy (weak), I3 = Ix (strong) — identical argument order to the
    # optimized global_K!(::Element{FixedFixed}) path above.
    local_stiffness_natural!(K, E, G, A, Iy, Ix, J, Ay, Az, L, aN, DN)
    # Symmetrize to remove any round-off asymmetry from the triple product.
    return 0.5 .* (K .+ transpose(K))
end

"""
    _condense_released_dofs(K, released)

Static condensation of the `released` local DOFs from a 12x12 stiffness matrix
`K`, returning a 12x12 matrix whose released rows/columns are zero and whose
retained block is `K_rr - K_rc K_cc^{-1} K_cr` (see derivation above). `K_cc`
must be non-singular (true for every release handled here).
"""
function _condense_released_dofs(K::Matrix{Float64}, released::Vector{Int})
    n = size(K, 1)
    retained = setdiff(1:n, released)
    Krr = K[retained, retained]
    Krc = K[retained, released]
    Kcc = K[released, released]
    # K_rc K_cc^{-1} K_cr, with K_cr = K_rc' by symmetry of K.
    Kred = Krr .- Krc * (Kcc \ transpose(Krc))
    out = zeros(n, n)
    out[retained, retained] .= 0.5 .* (Kred .+ transpose(Kred))
    # Snap catastrophic-cancellation round-off to exactly zero. When an entire
    # plane's rotations are released (pin-pin), the retained transverse-
    # translation stiffness is analytically zero (McGuire §7.4 Schur complement),
    # but the difference of two large near-equal terms leaves ~1e-18-relative
    # dust. The threshold is scaled to the matrix and sits many orders below any
    # physical stiffness coefficient, so no retained term is altered.
    tol = 1e-9 * maximum(abs, K)
    @inbounds for j in 1:n, i in 1:n
        if abs(out[i, j]) < tol
            out[i, j] = 0.0
        end
    end
    return out
end

"""
    k_fixedfixed(E, A, L, G, Ix, Iy, J; Ay = Inf, Az = Inf)

Full fixed-fixed 12x12 local stiffness matrix (no releases). Retains all six
end actions at both nodes. `Ay`/`Az` default to `Inf`, selecting the
Bernoulli-Euler formulation; finite values select Timoshenko shear flexibility.
"""
function k_fixedfixed(E, A, L, G, Ix, Iy, J; Ay = Inf, Az = Inf)
    return _k_fixedfixed_full(E, A, L, G, Ix, Iy, J, Ay, Az)
end

"""
    k_fixedfree(E, A, L, Ix, Iy)

Local stiffness for a member fixed at the start node i and fully rotationally
released at the end node j (release `:fixedfree`). Both bending planes become
fixed-pinned ("propped cantilever") members; releasing torsion at end j removes
torsional stiffness member-wide, so `G`/`J` are not required.

Built by condensing the released bending rotations theta_y,j (DOF 11) and
theta_z,j (DOF 12) from a fixed-fixed matrix assembled with G = J = 0 (torsion
DOFs 4 and 10 therefore carry zero rows/columns, matching the released torsion).
Kassimali "Structural Analysis" 5e, Table 8.2 (member with a hinge at one end):
transverse stiffness 3EI/L^3 and rotational stiffness 3EI/L at the fixed end.
"""
function k_fixedfree(E, A, L, Ix, Iy)
    Kff = _k_fixedfixed_full(E, A, L, 0.0, Ix, Iy, 0.0, Inf, Inf)
    return _condense_released_dofs(Kff, [11, 12])
end

"""
    k_freefixed(E, A, L, Ix, Iy)

Local stiffness for a member fully rotationally released at the start node i and
fixed at the end node j (release `:freefixed`) — the mirror image of
`k_fixedfree`. Bending planes become pinned-fixed members; torsion released at
start i removes torsional stiffness member-wide.

Built by condensing the released bending rotations theta_y,i (DOF 5) and
theta_z,i (DOF 6) from a fixed-fixed matrix assembled with G = J = 0.
Kassimali "Structural Analysis" 5e, Table 8.2 (member with a hinge at one end).
"""
function k_freefixed(E, A, L, Ix, Iy)
    Kff = _k_fixedfixed_full(E, A, L, 0.0, Ix, Iy, 0.0, Inf, Inf)
    return _condense_released_dofs(Kff, [5, 6])
end

"""
    k_freefree(E, A, L, Ix, Iy)

Local stiffness for a member with all rotational DOFs (both bending planes and
torsion) released at both ends (release `:freefree`). This is the classical
two-force / truss member: only axial stiffness EA/L survives on the axial DOFs
(1 and 7); every other row/column is zero.

Built by condensing both bending rotations of both planes — theta_y,i (DOF 5),
theta_z,i (DOF 6), theta_y,j (DOF 11), theta_z,j (DOF 12) — from a fixed-fixed
matrix assembled with G = J = 0 (so torsion DOFs 4 and 10 already carry zero
rows/columns, matching the released torsion at both ends). Condensing both
bending rotations of a plane drives that plane's transverse stiffness to exactly
zero (Kassimali "Structural Analysis" 5e — pin-pin member), leaving only the
axial two-force member.
"""
function k_freefree(E, A, L, Ix, Iy)
    Kff = _k_fixedfixed_full(E, A, L, 0.0, Ix, Iy, 0.0, Inf, Inf)
    return _condense_released_dofs(Kff, [5, 6, 11, 12])
end

"""
    k_joist(E, A, L, G, J, Ix, Iy)

Local stiffness for a joist / open-web-steel-joist idealization (release
`:joist`): bending is released at both ends in both planes, but axial (EA) and
torsion (GJ) are retained. This is the open-web-steel-joist idealization used
for OWSJ members.

Built by condensing the four released bending rotations — theta_y,i (DOF 5),
theta_z,i (DOF 6), theta_y,j (DOF 11), theta_z,j (DOF 12) — from a fixed-fixed
matrix assembled with the real G and J (torsion DOFs 4 and 10 stay in the
retained set). Static condensation drives both transverse (bending) planes to
exactly zero stiffness (Kassimali "Structural Analysis" 5e — pin-pin member), so
the result is the rank-2 matrix with axial EA/L on DOFs (1, 7) and St. Venant
torsion GJ/L on DOFs (4, 10); every other row/column is zero.
"""
function k_joist(E, A, L, G, J, Ix, Iy)
    Kff = _k_fixedfixed_full(E, A, L, G, Ix, Iy, J, Inf, Inf)
    return _condense_released_dofs(Kff, [5, 6, 11, 12])
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
    to_meters(element.length),
    to_meters_fourth(element.section.Ix),
    to_meters_fourth(element.section.Iy)
)

local_K(element::Element{Joist}) = k_joist(
    to_pascals(element.section.E),
    to_meters_squared(element.section.A),
    to_meters(element.length),
    to_pascals(element.section.G),
    to_meters_fourth(element.section.J),
    to_meters_fourth(element.section.Ix),
    to_meters_fourth(element.section.Iy)
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
