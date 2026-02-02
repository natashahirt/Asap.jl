#=
Geometric Stiffness for Frame Elements
======================================

Computes the geometric stiffness matrix (Kg) that accounts for the destabilizing
effect of axial loads on transverse stiffness.

Mathematical Background:
- Linear stiffness K captures material behavior
- Geometric stiffness Kg captures the P-Δ effect
- Total tangent stiffness: Kt = K + Kg
- Buckling: solve (K + λ·Kg)φ = 0

For a beam element with axial force P:
- Positive P (tension) increases transverse stiffness
- Negative P (compression) decreases transverse stiffness
- At critical load: det(K + Kg) = 0 → buckling

References:
- McGuire, Gallagher, Ziemian "Matrix Structural Analysis" Ch. 9
- Cook et al. "Concepts and Applications of FEA" Ch. 14
- Przemieniecki "Theory of Matrix Structural Analysis"
=#

using LinearAlgebra: norm
using SparseArrays: spzeros, SparseMatrixCSC

# =============================================================================
# Element-Level Geometric Stiffness
# =============================================================================

"""
    local_geometric_stiffness(L::Float64, P::Float64) -> Matrix{Float64}

Compute the 12×12 local geometric stiffness matrix for a frame element.

# Arguments
- `L::Float64`: Element length [m]
- `P::Float64`: Axial force [N] (positive = tension, negative = compression)

# Returns
- `Kg_local::Matrix{Float64}`: 12×12 geometric stiffness in local coordinates

# Local DOF Order
`[u₁, v₁, w₁, θx₁, θy₁, θz₁, u₂, v₂, w₂, θx₂, θy₂, θz₂]`

# Notes
The geometric stiffness is proportional to P/L and affects only the
transverse DOFs (not axial or torsional).

Based on the consistent geometric stiffness formulation from
Przemieniecki (1968) and Cook et al.
"""
function local_geometric_stiffness(L::Float64, P::Float64)
    # Geometric stiffness only affects transverse behavior
    # Standard formulation from Przemieniecki / Cook et al.
    
    Kg = zeros(12, 12)
    
    # Factor
    c = P / L
    
    # Coefficients for consistent geometric stiffness
    # These come from the virtual work of the axial force through
    # the second-order transverse displacement terms
    
    a = 6.0 / 5.0        # Diagonal terms for translations
    b = L / 10.0         # Coupling terms for rotation-translation  
    d = 2.0 * L^2 / 15.0 # Diagonal terms for rotations
    e = -L^2 / 30.0      # Off-diagonal rotation terms
    
    # Local DOF indices:
    # 1=u1, 2=v1, 3=w1, 4=θx1, 5=θy1, 6=θz1
    # 7=u2, 8=v2, 9=w2, 10=θx2, 11=θy2, 12=θz2
    
    # v-direction (local y, bending about z)
    # v1-v1
    Kg[2, 2] = a * c
    # v1-θz1
    Kg[2, 6] = b * c
    Kg[6, 2] = b * c
    # v1-v2
    Kg[2, 8] = -a * c
    Kg[8, 2] = -a * c
    # v1-θz2
    Kg[2, 12] = b * c
    Kg[12, 2] = b * c
    # θz1-θz1
    Kg[6, 6] = d * c
    # θz1-v2
    Kg[6, 8] = -b * c
    Kg[8, 6] = -b * c
    # θz1-θz2
    Kg[6, 12] = e * c
    Kg[12, 6] = e * c
    # v2-v2
    Kg[8, 8] = a * c
    # v2-θz2
    Kg[8, 12] = -b * c
    Kg[12, 8] = -b * c
    # θz2-θz2
    Kg[12, 12] = d * c
    
    # w-direction (local z, bending about y)
    # Same structure but with opposite sign for coupling terms
    # w1-w1
    Kg[3, 3] = a * c
    # w1-θy1 (note: opposite sign due to rotation convention)
    Kg[3, 5] = -b * c
    Kg[5, 3] = -b * c
    # w1-w2
    Kg[3, 9] = -a * c
    Kg[9, 3] = -a * c
    # w1-θy2
    Kg[3, 11] = -b * c
    Kg[11, 3] = -b * c
    # θy1-θy1
    Kg[5, 5] = d * c
    # θy1-w2
    Kg[5, 9] = b * c
    Kg[9, 5] = b * c
    # θy1-θy2
    Kg[5, 11] = e * c
    Kg[11, 5] = e * c
    # w2-w2
    Kg[9, 9] = a * c
    # w2-θy2
    Kg[9, 11] = b * c
    Kg[11, 9] = b * c
    # θy2-θy2
    Kg[11, 11] = d * c
    
    return Kg
end

"""
    geometric_stiffness(element::Element) -> Matrix{Float64}

Compute the 12×12 global geometric stiffness matrix for a frame element
using its current internal axial force.

# Returns
- `Kg::Matrix{Float64}`: 12×12 geometric stiffness in global coordinates

# Notes
- Uses the axial force from the most recent solve
- Call `solve!(model)` first to get the internal forces
- The element's transformation matrix R is used for coordinate transformation
"""
function geometric_stiffness(element::FrameElement)
    # Get axial force from element (in local coordinates)
    # element.forces[1] is axial force at node 1 (tension positive)
    P = -element.forces[1]  # Use average, negative because internal compression
    
    # Get length
    L = ustrip(u"m", element.length)
    
    # Compute local geometric stiffness
    Kg_local = local_geometric_stiffness(L, P)
    
    # Transform to global coordinates: Kg_global = R' * Kg_local * R
    R = element.R
    Kg_global = R' * Kg_local * R
    
    return Kg_global
end

"""
    geometric_stiffness(element::Element, P::Float64) -> Matrix{Float64}

Compute the 12×12 global geometric stiffness matrix for a frame element
using a specified axial force.

# Arguments
- `element::Element`: The frame element
- `P::Float64`: Axial force in Newtons (positive = tension)
"""
function geometric_stiffness(element::FrameElement, P::Float64)
    L = ustrip(u"m", element.length)
    Kg_local = local_geometric_stiffness(L, P)
    R = element.R
    return R' * Kg_local * R
end

# =============================================================================
# Global Assembly
# =============================================================================

"""
    assemble_geometric_stiffness(model::AbstractModel) -> SparseMatrixCSC

Assemble the global geometric stiffness matrix from current internal forces.

# Prerequisites
- Model must be processed (`process!(model)`)
- Static analysis must be run (`solve!(model)`) to compute internal forces

# Returns
- `Kg::SparseMatrixCSC{Float64}`: Global geometric stiffness matrix

# Example
```julia
model = Model(nodes, elements, loads)
solve!(model)  # Get internal forces first
Kg = assemble_geometric_stiffness(model)
```
"""
function assemble_geometric_stiffness(model::FrameModel)
    !model.processed && error("Model must be processed. Call process!(model) first.")
    
    n_dof = model.nDOFs
    Kg = spzeros(Float64, n_dof, n_dof)
    
    for element in model.elements
        # Skip non-frame elements or elements without forces
        if !(element isa FrameElement)
            continue
        end
        
        # Compute element geometric stiffness
        Kg_elem = geometric_stiffness(element)
        
        # Assemble into global matrix
        gid = element.globalID
        for i in 1:12
            for j in 1:12
                Kg[gid[i], gid[j]] += Kg_elem[i, j]
            end
        end
    end
    
    return Kg
end

function assemble_geometric_stiffness(model::Model)
    !model.processed && error("Model must be processed. Call process!(model) first.")
    
    n_dof = model.nDOFs
    Kg = spzeros(Float64, n_dof, n_dof)
    
    # Frame elements
    for element in model.frame_elements
        if !(element isa FrameElement)
            continue
        end
        
        Kg_elem = geometric_stiffness(element)
        gid = element.globalID
        for i in 1:12
            for j in 1:12
                Kg[gid[i], gid[j]] += Kg_elem[i, j]
            end
        end
    end
    
    # Shell elements
    for elem in model.shell_elements
        Kg_elem = geometric_stiffness(elem, model)
        gid = elem.globalID
        for i in 1:18
            for j in 1:18
                Kg[gid[i], gid[j]] += Kg_elem[i, j]
            end
        end
    end
    
    return Kg
end

# =============================================================================
# Shell Element Geometric Stiffness
# =============================================================================

"""
    local_geometric_stiffness(elem::ShellTri3, σ_membrane::Vector{Float64}) -> Matrix{Float64}

Compute the 18×18 local geometric stiffness matrix for a ShellTri3 element.

# Mathematical Background
The geometric stiffness for a flat shell under membrane forces (Nxx, Nyy, Nxy)
captures the destabilizing effect of in-plane loads on out-of-plane stability.

The geometric potential energy is:
    Πg = (1/2) ∫∫ [Nxx(∂w/∂x)² + Nyy(∂w/∂y)² + 2Nxy(∂w/∂x)(∂w/∂y)] dA

This leads to:
    Kg = Ae · Bg' · σ · Bg

where:
- Bg relates nodal DOFs to ∂w/∂x, ∂w/∂y (using shape function gradients)
- σ = [Nxx Nxy; Nxy Nyy] is the membrane force tensor
- Ae is the element area

# Arguments
- `elem::ShellTri3`: Processed shell element (must have ecoords_e computed)
- `σ_membrane::Vector{Float64}`: [Nxx, Nyy, Nxy] membrane forces [N/m]

# Returns
- `Kg::Matrix{Float64}`: 18×18 geometric stiffness in local coordinates

# References
- Cook et al. "Concepts and Applications of FEA" Ch. 14
- Przemieniecki "Theory of Matrix Structural Analysis"
- FEniCSx-Shells nonlinear Naghdi formulation
"""
function local_geometric_stiffness(elem::ShellTri3, σ_membrane::Vector{Float64})
    Nxx, Nyy, Nxy = σ_membrane[1], σ_membrane[2], σ_membrane[3]
    Ae = elem.area
    
    # Compute shape function gradients in local coordinates
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    # Build Bg matrix: relates nodal w displacements to [∂w/∂x; ∂w/∂y]
    # For CST-type element: w = Σ Ni * wi, so ∂w/∂x = Σ (∂Ni/∂x) * wi
    # Bg is 2×3 (2 gradients, 3 nodes with w DOFs)
    Bg = zeros(2, 3)
    for i in 1:3
        Bg[1, i] = gradN_e[i, 1]  # ∂Ni/∂x
        Bg[2, i] = gradN_e[i, 2]  # ∂Ni/∂y
    end
    
    # Membrane stress tensor (2×2)
    σ_mat = [Nxx Nxy; Nxy Nyy]
    
    # Geometric stiffness for w DOFs: Kg_ww = Ae · Bg' · σ · Bg (3×3)
    Kg_ww = Ae * (Bg' * σ_mat * Bg)
    
    # Expand to full 18×18 element matrix
    # DOF ordering per node: [u, v, w, θx, θy, θz]
    # Node 1: DOFs 1-6, Node 2: DOFs 7-12, Node 3: DOFs 13-18
    # w DOFs are at indices 3, 9, 15
    Kg = zeros(18, 18)
    w_dofs = [3, 9, 15]
    
    for i in 1:3
        for j in 1:3
            Kg[w_dofs[i], w_dofs[j]] = Kg_ww[i, j]
        end
    end
    
    return Kg
end

"""
    geometric_stiffness(elem::ShellTri3, model::AbstractModel) -> Matrix{Float64}

Compute the 18×18 global geometric stiffness matrix for a shell element
using its current membrane forces from a solved model.

# Arguments
- `elem::ShellTri3`: The shell element
- `model::AbstractModel`: Solved model with displacement vector

# Returns
- `Kg_global::Matrix{Float64}`: 18×18 geometric stiffness in global coordinates
"""
function geometric_stiffness(elem::ShellTri3, model::AbstractModel)
    # Get membrane forces from solved model
    sif = ShellInternalForces(elem, model.u)
    σ_membrane = [sif.Nxx, sif.Nyy, sif.Nxy]
    
    # Compute local geometric stiffness
    Kg_local = local_geometric_stiffness(elem, σ_membrane)
    
    # Transform to global coordinates: Kg_global = R' · Kg_local · R
    Kg_global = elem.R' * Kg_local * elem.R
    
    return Kg_global
end

"""
    geometric_stiffness(elem::ShellTri3, σ_membrane::Vector{Float64}) -> Matrix{Float64}

Compute geometric stiffness with specified membrane forces.

# Arguments
- `elem::ShellTri3`: Processed shell element
- `σ_membrane::Vector{Float64}`: [Nxx, Nyy, Nxy] membrane forces [N/m]
"""
function geometric_stiffness(elem::ShellTri3, σ_membrane::Vector{Float64})
    Kg_local = local_geometric_stiffness(elem, σ_membrane)
    return elem.R' * Kg_local * elem.R
end

# Support for CompositeShellTri3
function local_geometric_stiffness(elem::CompositeShellTri3, σ_membrane::Vector{Float64})
    # Same formulation - membrane forces affect out-of-plane stability identically
    Nxx, Nyy, Nxy = σ_membrane[1], σ_membrane[2], σ_membrane[3]
    Ae = elem.area
    
    gradN_e = zeros(3, 2)
    gradN_e, _ = _compute_gradN_Ae!(gradN_e, elem.ecoords_e)
    
    Bg = zeros(2, 3)
    for i in 1:3
        Bg[1, i] = gradN_e[i, 1]
        Bg[2, i] = gradN_e[i, 2]
    end
    
    σ_mat = [Nxx Nxy; Nxy Nyy]
    Kg_ww = Ae * (Bg' * σ_mat * Bg)
    
    Kg = zeros(18, 18)
    w_dofs = [3, 9, 15]
    for i in 1:3
        for j in 1:3
            Kg[w_dofs[i], w_dofs[j]] = Kg_ww[i, j]
        end
    end
    
    return Kg
end

function geometric_stiffness(elem::CompositeShellTri3, model::AbstractModel)
    sif = ShellInternalForces(elem, model.u)
    σ_membrane = [sif.Nxx, sif.Nyy, sif.Nxy]
    Kg_local = local_geometric_stiffness(elem, σ_membrane)
    return elem.R' * Kg_local * elem.R
end

# =============================================================================
# Assembly for Shell Models
# =============================================================================

"""
    assemble_geometric_stiffness(model::ShellModel) -> SparseMatrixCSC

Assemble global geometric stiffness matrix for shell-only models.
"""
function assemble_geometric_stiffness(model::ShellModel)
    !model.processed && error("Model must be processed. Call process!(model) first.")
    
    n_dof = model.nDOFs
    Kg = spzeros(Float64, n_dof, n_dof)
    
    for elem in model.elements
        Kg_elem = geometric_stiffness(elem, model)
        gid = elem.globalID
        for i in 1:18
            for j in 1:18
                Kg[gid[i], gid[j]] += Kg_elem[i, j]
            end
        end
    end
    
    return Kg
end

"""
    assemble_geometric_stiffness(model::ShellModel, σ_uniform::Vector{Float64}) -> SparseMatrixCSC

Assemble global geometric stiffness matrix with specified uniform membrane forces.

This is useful for buckling analysis when you want to prescribe a known stress state
(e.g., for validation against analytical solutions) rather than computing it from
a static solve.

# Arguments
- `model::ShellModel`: Processed shell model
- `σ_uniform::Vector{Float64}`: [Nxx, Nyy, Nxy] uniform membrane forces [N/m]

# Example
```julia
# Uniform uniaxial compression of 1000 N/m
Kg = assemble_geometric_stiffness(model, [-1000.0, 0.0, 0.0])
```
"""
function assemble_geometric_stiffness(model::ShellModel, σ_uniform::Vector{Float64})
    !model.processed && error("Model must be processed. Call process!(model) first.")
    
    n_dof = model.nDOFs
    Kg = spzeros(Float64, n_dof, n_dof)
    
    for elem in model.elements
        Kg_elem = geometric_stiffness(elem, σ_uniform)
        gid = elem.globalID
        for i in 1:18
            for j in 1:18
                Kg[gid[i], gid[j]] += Kg_elem[i, j]
            end
        end
    end
    
    return Kg
end

# TrussModel - simplified for axial-only elements
function assemble_geometric_stiffness(model::TrussModel)
    !model.processed && error("Model must be processed. Call process!(model) first.")
    
    n_dof = model.nDOFs
    Kg = spzeros(Float64, n_dof, n_dof)
    
    for element in model.elements
        # For truss elements, geometric stiffness is simpler
        # Kg = (P/L) * [c² cs; cs s²] where c,s are direction cosines
        P = -element.forces[1]  # Axial force (compression negative)
        L = ustrip(u"m", element.length)
        
        # Direction cosines
        n1_pos = element.nodeStart.position
        n2_pos = element.nodeEnd.position
        dx = [ustrip(u"m", n2_pos[i] - n1_pos[i]) for i in 1:3]
        
        c = P / L
        
        # Build 6×6 geometric stiffness for 3D truss
        # Affects transverse DOFs relative to element axis
        Kg_elem = zeros(6, 6)
        
        # For truss, geometric stiffness in each perpendicular direction
        for i in 1:3
            for j in 1:3
                if i == j
                    # Diagonal: (1 - n_i²) where n is unit vector along element
                    ni = dx[i] / L
                    val = c * (1 - ni^2)
                else
                    # Off-diagonal: -n_i * n_j
                    ni = dx[i] / L
                    nj = dx[j] / L
                    val = -c * ni * nj
                end
                Kg_elem[i, i] += val / 2  # Factor of 2 for assembly
                Kg_elem[i+3, i+3] += val / 2
                Kg_elem[i, i+3] -= val / 2
                Kg_elem[i+3, i] -= val / 2
            end
        end
        
        # Assemble
        gid = element.globalID
        for i in 1:6
            for j in 1:6
                Kg[gid[i], gid[j]] += Kg_elem[i, j]
            end
        end
    end
    
    return Kg
end
