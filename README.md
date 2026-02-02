[![DOI](https://zenodo.org/badge/426740094.svg)](https://zenodo.org/doi/10.5281/zenodo.10581559)

![](READMEassets/forces-axo.png)

# Asap.jl

Asap is...

- Another Structural Analysis Package
- results As Soon As Possible
- Analysis of Structures Avec Programming
- A Simple Analysis Please

Designed first-and-foremost for information-rich data structures and ease of querying, but always with performance in mind.

## About This Fork

This is an extended fork of the [original Asap.jl](https://github.com/keithjlee/Asap) by Keith Janghyun Lee. It adds significant new capabilities while maintaining API compatibility:

**New in this fork:**

- **Shell elements** - Triangular shell elements with membrane + bending (based on FinEtools)
- **Composite shells** - Laminated composite materials with ply-level stress recovery
- **Grounded springs** - Elastic foundation/support modeling
- **Modal analysis** - Natural frequencies and mode shapes via `modal!`
- **Mixed models** - Combined frame + shell structures in one model
- **Tributary areas** - Straight skeleton and Voronoi algorithms for load distribution
- **Area loads** - Unified surface pressure API for shells
- **Force Density Method** - Form-finding for cable/membrane structures

### Installation

Since this is a local fork, add it as a development dependency:

```julia
using Pkg
Pkg.develop(path="path/to/Asap")
```

Or in package mode:

```julia
pkg> dev path/to/Asap
```

---

# Quick Start

Here's a complete example: a 2-story, 2×2 bay building with concrete columns, steel beams, shell floor slabs, foundation springs, and modal analysis.

```julia
using Asap
using Unitful

# =============================================================================
# Materials
# =============================================================================

# Frame sections
concrete_col = Section(
    0.09u"m^2",      # 300×300mm column
    30u"GPa",        # E (concrete)
    12.5u"GPa",      # G
    6.75e-4u"m^4",   # Ix
    6.75e-4u"m^4",   # Iy
    1.0e-3u"m^4"     # J
)

steel_beam = Section(
    7.6e-3u"m^2",    # W310×33
    200u"GPa",       # E (steel)
    77u"GPa",        # G
    6.5e-5u"m^4",    # Ix
    1.4e-5u"m^4",    # Iy
    2.0e-7u"m^4"     # J
)

# Shell material
slab_material = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)

# =============================================================================
# Geometry: 2×2 bays, 6m × 6m each, 2 stories @ 3.5m
# =============================================================================

Lx, Ly = 6.0u"m", 6.0u"m"    # Bay dimensions
H = 3.5u"m"                   # Story height
nx, ny = 3, 3                 # 3×3 grid of columns

# Create column nodes (ground level fixed, upper levels free)
nodes = Node[]
for k in 0:2  # 3 levels: ground, level 1, level 2
    for j in 0:ny-1
        for i in 0:nx-1
            pos = [i*Lx, j*Ly, k*H]
            fixity = k == 0 ? :fixed : :free
            push!(nodes, Node(pos, fixity))
        end
    end
end

# Helper to get node index
node_idx(i, j, k) = k * nx * ny + j * nx + i + 1

# =============================================================================
# Frame Elements: Columns and Beams
# =============================================================================

elements = Element[]

# Columns (vertical)
for k in 0:1  # Between levels
    for j in 0:ny-1
        for i in 0:nx-1
            n_bot = nodes[node_idx(i, j, k)]
            n_top = nodes[node_idx(i, j, k+1)]
            push!(elements, Element(n_bot, n_top, concrete_col, :column))
        end
    end
end

# Beams (X-direction)
for k in 1:2  # At each floor level
    for j in 0:ny-1
        for i in 0:nx-2
            n1 = nodes[node_idx(i, j, k)]
            n2 = nodes[node_idx(i+1, j, k)]
            push!(elements, Element(n1, n2, steel_beam, :beam_x))
        end
    end
end

# Beams (Y-direction)
for k in 1:2
    for j in 0:ny-2
        for i in 0:nx-1
            n1 = nodes[node_idx(i, j, k)]
            n2 = nodes[node_idx(i, j+1, k)]
            push!(elements, Element(n1, n2, steel_beam, :beam_y))
        end
    end
end

# =============================================================================
# Shell Floor Slabs
# =============================================================================

# Define shell section (thickness + material)
slab_section = ShellSection(0.15u"m", slab_material)

shell_elements = ShellTri3[]

for k in 1:2  # Each floor level
    for j in 0:ny-2
        for i in 0:nx-2
            # Get bay corners (counter-clockwise)
            corners = (
                nodes[node_idx(i, j, k)],
                nodes[node_idx(i+1, j, k)],
                nodes[node_idx(i+1, j+1, k)],
                nodes[node_idx(i, j+1, k)]
            )
            level_id = k == 1 ? :floor1 : :floor2
            
            # Shell() auto-triangulates any polygon!
            append!(shell_elements, Shell(corners, slab_section; id=level_id))
        end
    end
end

# Get all nodes from shells (includes corners + interior mesh nodes)
shell_nodes = get_nodes(shell_elements)

# =============================================================================
# Foundation Springs (at base nodes)
# =============================================================================

springs = Spring[]
k_vertical = 50000u"kN/m"    # Vertical soil spring
k_horizontal = 10000u"kN/m"  # Horizontal soil spring

for j in 0:ny-1
    for i in 0:nx-1
        base_node = nodes[node_idx(i, j, 0)]
        push!(springs, Spring(base_node; 
            kx=k_horizontal, ky=k_horizontal, kz=k_vertical))
    end
end

# =============================================================================
# Loads
# =============================================================================

# Floor live load as area load
floor_pressure = 2.5u"kPa"
loads = AbstractLoad[]

for shell in shell_elements
    push!(loads, AreaLoad([shell], floor_pressure))
end

# =============================================================================
# Model and Analysis
# =============================================================================

# Create unified model (frames + shells)
model = Model(nodes, elements, shell_elements, loads)

# Solve static analysis
solve!(model)

# Add foundation springs and re-solve
process!(model)
add_springs!(model, springs)
solve!(model)

# Modal analysis
result = modal!(model; n=10)
print_modal_summary(result)

# Access results
println("\\nFirst 3 natural frequencies: ", round.(result.frequencies[1:3], digits=2), " Hz")
println("Base reaction at corner: ", nodes[1].reaction)
```

---

# Units

Asap uses [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) for type-safe unit handling. All quantities convert to SI internally (meters, Newtons, Pascals).

```julia
using Unitful

# Any compatible units work
node = Node([10.0u"ft", 20.0u"ft", 0.0u"inch"], :fixed)
load = NodeForce(node, [0.0u"kip", -50.0u"kip", 0.0u"kip"])

# Results return in SI; convert as needed
d_mm = uconvert(u"mm", node.displacement[2])
R_kip = uconvert(u"kip", node.reaction[2])
```

---

# Core Types

## Nodes

Nodes define spatial positions and boundary conditions.

### `Node` (6 DOF)

For frame and shell models: 3 translations + 3 rotations.

```julia
# Explicit DOFs: [Tx, Ty, Tz, Rx, Ry, Rz]
node = Node([0.0u"m", 0.0u"m", 0.0u"m"], [false, false, false, true, true, true])

# Using fixity symbols
node = Node([0.0u"m", 5.0u"m", 0.0u"m"], :fixed)   # All DOFs restrained
node = Node([0.0u"m", 5.0u"m", 0.0u"m"], :pinned)  # Translations fixed, rotations free
node = Node([0.0u"m", 5.0u"m", 0.0u"m"], :free)    # All DOFs free

# With identifier (for filtering)
node = Node([0.0u"m", 5.0u"m", 0.0u"m"], :free, :floor_node)
```

Available fixity symbols: `:fixed`, `:free`, `:pinned`, `:xfixed`, `:yfixed`, `:zfixed`, `:xfree`, `:yfree`, `:zfree`

### `TrussNode` (3 DOF)

For pure truss models: translations only (optimized, smaller matrices).

```julia
node = TrussNode([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
```

> **Note:** Use `TrussNode` with `TrussElement` and `TrussModel` for pure truss structures. For truss-like elements in mixed models, use regular `Node` with `Element(...; release=:freefree)`.

## Springs

Grounded springs connect nodes to "ground" with elastic stiffness. Useful for foundation modeling, elastic supports, and isolators.

```julia
# Uniform translational stiffness
spring = Spring(node, 1000u"kN/m")

# Per-DOF specification
spring = Spring(node; kx=100u"kN/m", ky=100u"kN/m", kz=500u"kN/m")

# With rotational stiffness
spring = Spring(node; kz=500u"kN/m", krz=1000u"kN*m/rad")

# Add springs after processing
process!(model)
add_springs!(model, springs)
solve!(model)
```

## Sections & Materials

### `Section` (Frame Elements)

Defines cross-section properties for beams and columns.

```julia
section = Section(A, E, G, Ix, Iy, J)

# Example: W310×33 steel beam
w310 = Section(
    4.18e-3u"m^2",   # A: Area
    200u"GPa",       # E: Young's modulus
    77u"GPa",        # G: Shear modulus
    6.5e-5u"m^4",    # Ix: Strong axis moment of inertia
    1.4e-5u"m^4",    # Iy: Weak axis moment of inertia
    2.0e-7u"m^4"     # J: Torsional constant
)
```

### `TrussSection` (Truss Elements)

Simplified section for axial-only members.

```julia
section = TrussSection(A, E)
section = TrussSection(0.01u"m^2", 200u"GPa")
```

### `ShellMaterial` (Shell Elements)

Material definition for shells. Iteration-friendly for parametric studies.

```julia
concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
steel = ShellMaterial(E=200u"GPa", ν=0.3, ρ=7850u"kg/m^3", name=:steel)

# Built-in presets
Concrete_Shell  # 30 GPa, ν=0.2, ρ=2400 kg/m³
Steel_Shell     # 200 GPa, ν=0.3, ρ=7850 kg/m³

# Easy iteration
for mat in [concrete, steel, aluminum]
    shells = [ShellTri3(nodes, 0.15u"m", mat) for nodes in mesh]
    model = ShellModel(all_nodes, shells, loads)
    solve!(model)
end
```

### Composite Materials

For laminated shells, define plies and stacking sequences:

```julia
# Define a ply
carbon = Ply("T300/5208", E1=181u"GPa", E2=10.3u"GPa", G12=7.17u"GPa", 
             ν12=0.28, ρ=1600u"kg/m^3", t=0.125u"mm")

# Create laminate [0/90/±45]s
laminate = symmetric_laminate([
    (carbon, 0.0),
    (carbon, 90.0),
    (carbon, 45.0),
    (carbon, -45.0)
])

# Use in composite shell
shell = CompositeShellTri3((n1, n2, n3), laminate)
```

## Elements

### `Element` (Frame)

6-DOF frame element for beams and columns.

```julia
element = Element(node_start, node_end, section)
element = Element(node_start, node_end, section, :beam_id)

# With end releases (hinges)
element = Element(n1, n2, section; release=:freefree)  # Pin-pin (truss behavior)
element = Element(n1, n2, section; release=:fixedfree) # Fixed-pin
element = Element(n1, n2, section; release=:joist)     # Joist connection
```

Release options: `:fixedfixed` (default), `:fixedfree`, `:freefixed`, `:freefree`, `:joist`

### `TrussElement`

Axial-only element for truss structures. Use with `TrussNode` and `TrussModel`.

```julia
element = TrussElement(truss_node1, truss_node2, truss_section)
```

### `ShellTri3`

3-node triangular shell element with membrane, bending, and shear stiffness.

```julia
# With explicit properties
shell = ShellTri3((n1, n2, n3), 0.2u"m", 30u"GPa", 0.2; ρ=2400.0)

# With ShellMaterial (preferred for iteration)
shell = ShellTri3((n1, n2, n3), 0.2u"m", concrete_material)
```

### `CompositeShellTri3`

Laminated composite shell for layered materials.

```julia
shell = CompositeShellTri3((n1, n2, n3), laminate; id=:wing_skin)
```

### `ShellSection`

Combines thickness with material properties. Use this instead of passing thickness and material separately - it's cleaner and lets you bring your own material definitions.

```julia
# Direct specification
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete)

# From ShellMaterial
concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3")
section = ShellSection(0.15u"m", concrete)

# Built-in presets
Concrete_ShellSection_150mm  # 150mm concrete slab
Concrete_ShellSection_200mm  # 200mm concrete slab
```

### `Shell()` - Create Shell Elements

Automatically triangulate any polygon into `ShellTri3` elements. Works with triangles, quads, or any convex/concave polygon.

```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")

# Simple - uses default n=4 refinement (~5% accuracy)
shells = Shell((n1, n2, n3, n4), section)

# Custom refinement (n=6 gives ~1% accuracy)
shells = Shell((n1, n2, n3, n4), 6, section)

# Custom element ID
shells = Shell(corners, section; id=:floor_slab)
```

Use `get_nodes(shells)` to extract all nodes (corners + interior) for the model:
```julia
shells = Shell(corners, section)
model = ShellModel(get_nodes(shells), shells, loads)
```

#### Interior Supports (Continuous Slabs)

For slabs spanning multiple bays with interior beams, use `interior_supports` to ensure correct moment distribution:

```julia
# Multi-bay slab with interior support beam
shells = Shell((n1, n2, n3, n4), section;
    interior_supports = [interior_beam],      # Elements or (Node, Node) pairs
    edge_support_type = :pinned,              # Fixity for perimeter nodes
    interior_support_type = :pinned           # Fixity for interior support nodes
)
```

This automatically:
1. Adds mesh nodes along the interior support line
2. Applies the specified fixity to those nodes
3. Ensures proper moment continuity over interior supports

**Support types:**
| Symbol | Description | Fixed DOFs |
|--------|-------------|------------|
| `:pinned` | Simply-supported | x, y, z translations |
| `:zfixed` | Vertical support only | z translation |
| `:fixed` | Fully clamped | All 6 DOFs |
| `:free` | No constraint | None |
| `Vector{Bool}` | Custom | `[x, y, z, θx, θy, θz]` |

**Example: Two-span continuous slab**
```julia
# Define geometry (8m x 4m slab, interior beam at x=4m)
n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
n2 = Node([8.0u"m", 0.0u"m", 0.0u"m"], :pinned)
n3 = Node([8.0u"m", 4.0u"m", 0.0u"m"], :pinned)
n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)

# Interior beam at midspan
beam_n1 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
beam_n2 = Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
interior_beam = Element(beam_n1, beam_n2, beam_section)

# Create continuous slab
sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
shells = Shell((n1, n2, n3, n4), 6, sec;
    interior_supports = [interior_beam]
)

# Analysis
model = ShellModel(get_nodes(shells), shells, [AreaLoad(shells, 5u"kPa")])
solve!(model)
```

The interior support creates pinned nodes along the beam line, producing correct two-span moment distribution (hogging moment over support, reduced span moments).

### `mesh()` - Just Triangulation

If you want to control meshing separately:
```julia
tri = mesh((n1, n2, n3, n4), 6)  # Returns DelaunayTriangulation object
shells = Shell(tri, section)     # Create shells from triangulation
```

---

# Loads

All load values require Unitful units.

## `NodeForce`

Point force at a node.

```julia
load = NodeForce(node, [0.0u"kN", 0.0u"kN", -50.0u"kN"])
```

## `NodeMoment`

Point moment at a node.

```julia
load = NodeMoment(node, [0.0u"kN*m", 100.0u"kN*m", 0.0u"kN*m"])
```

## `LineLoad`

Distributed load along a frame element.

```julia
load = LineLoad(element, [0.0u"kN/m", 0.0u"kN/m", -10.0u"kN/m"])
```

## `PointLoad`

Concentrated load at a position along an element.

```julia
load = PointLoad(element, 0.5, [0.0u"kN", -20.0u"kN", 0.0u"kN"])  # At midspan
```

## `AreaLoad`

Surface pressure on shell elements. The unified API for floor/slab loading.

```julia
# FEM approach: forces at shell nodes (default)
load = AreaLoad(shells, 5.0u"kPa")

# Tributary distribution to edge beams
load = AreaLoad(shells, 5.0u"kPa"; distribute_to=edge_beams)

# Tributary with interior beams (interior beams get load from both sides)
load = AreaLoad(shells, 5.0u"kPa"; 
    distribute_to=vcat(edge_beams, interior_beam),
    interior_beams=[interior_beam]
)

# One-way slab (tributary distribution along X-axis)
load = AreaLoad(shells, 5.0u"kPa"; distribute_to=beams, axis=(1.0, 0.0))

# Custom load direction (default is -Z gravity)
load = AreaLoad(shells, 5.0u"kPa"; direction=(0.0, 0.0, -1.0))
```

## `SelfWeight`

Automatic self-weight from shell properties (`p = ρ × t × g`).

```julia
# Shell self-weight (requires ρ in ShellSection)
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
shells = Shell(corners, section)

sw = SelfWeight(shells)  # Computes: 2400 × 0.15 × 9.81 = 3531 Pa

# Custom gravity
sw = SelfWeight(shells; g=10u"m/s^2")
```

---

# Models

Asap provides specialized model types that auto-dispatch to optimized solvers.

## `Model` (Unified)

For mixed frame + shell structures. The most general model type.

```julia
model = Model(nodes, frame_elements, shell_elements, loads)

# Frame-only convenience (creates Model with empty shells)
model = Model(nodes, beams, loads)
```

## `FrameModel`

Frame elements only. Equivalent to `Model` without shells.

```julia
model = FrameModel(nodes, elements, loads)
```

## `ShellModel`

Shell elements only.

```julia
model = ShellModel(nodes, shell_elements, loads)
```

## `TrussModel`

Optimized for pure truss structures (3 DOF per node).

```julia
model = TrussModel(truss_nodes, truss_elements, node_forces)
```

## Diaphragm Behavior

Shell elements inherently provide **diaphragm action** through their membrane (in-plane) stiffness. When shell elements share nodes with frame elements, they automatically:

- Transfer lateral loads to the frame
- Distribute forces between columns
- Provide in-plane rigidity

**You don't need to add explicit diaphragm constraints**—the shell membrane stiffness handles it naturally through FEA.

```julia
# Concrete slab at floor level (E=30GPa, t=150mm)
# Provides very stiff diaphragm through membrane action
section = ShellSection(0.15u"m", 30u"GPa", 0.2)
floor_shells = Shell((col1, col2, col3, col4), section)

# Share nodes with columns → automatic load transfer
model = Model(all_nodes, columns, floor_shells, lateral_loads)
```

For **rigid diaphragm** assumptions (infinite in-plane stiffness), the shell stiffness is typically sufficient for concrete floors. For explicit rigid diaphragms without mass, use the `Diaphragm` type:

### `Diaphragm()` - Massless Rigid Diaphragm

Create shell elements that are very stiff in-plane but contribute no mass:

```julia
# Basic usage - connects frame nodes at floor level
diaphragms = Diaphragm((col1_top, col2_top, col3_top, col4_top))

# With custom stiffness
soft_diaphragm = DiaphragmSection(E=200u"GPa")  # Steel-like
diaphragms = Diaphragm(corners, soft_diaphragm)
```

**Properties:**
| Property | Default | Purpose |
|----------|---------|---------|
| `E` | 1 TPa | Very high → rigid in-plane |
| `ρ` | 0 | No mass → doesn't affect dynamics |
| `thickness` | 10mm | Small → minimal out-of-plane stiffness |

**Use `Diaphragm` when:**
- Lateral analysis of frame buildings
- You need rigid in-plane behavior but don't care about slab moments
- Mass should come only from frame elements (explicit mass modeling)

**Use `Shell` when:**
- You need slab bending moments for design
- Mass distribution matters (seismic/modal analysis)
- Accurate slab deflections needed

```julia
# Diaphragm for lateral load distribution
n1, n2, n3, n4 = [col.top_node for col in floor_columns]
diaphragms = Diaphragm((n1, n2, n3, n4))

model = Model(all_nodes, columns, diaphragms, lateral_loads)
solve!(model)  # Lateral loads distributed via rigid diaphragm
```

---

# Analysis

## Static Analysis

```julia
solve!(model)  # Process and solve in one call

# Access results
node.displacement  # [Tx, Ty, Tz, Rx, Ry, Rz] with units
node.reaction      # [Fx, Fy, Fz, Mx, My, Mz] with units
element.forces     # Local end forces

# Re-solve with new loads
solve!(model, new_loads)

# After modifying geometry
solve!(model; reprocess=true)
```

## Modal Analysis

Compute natural frequencies and mode shapes.

```julia
# Primary API (auto-processes if needed)
result = modal!(model)           # 6 modes (default)
result = modal!(model; n=10)     # 10 modes

# Access results
result.frequencies   # Natural frequencies [Hz]
result.periods       # Periods [s]
result.mode_shapes   # Mode shape matrix (columns = modes)
result.omegas        # Angular frequencies [rad/s]

# Pretty print summary
print_modal_summary(result)

# Mass matrix options
result = modal!(model; n=6, mass_type=MASS_LUMPED)
```

Mass types: `MASS_CONSISTENT` (default), `MASS_LUMPED`, `MASS_CONSISTENT_NO_ROTATION`, `MASS_LUMPED_NO_ROTATION`

---

# Post-Processing

Originally AsapToolkit features now native in Asap.

## Frame Internal Forces

Extract forces and moments along beam length.

```julia
forces = ElementInternalForces(element, model)

# At any position along element (0 to 1)
Vy = forces.Vy(0.5)   # Shear at midspan
Mz = forces.Mz(0.5)   # Moment at midspan

# Envelopes for multiple load cases
envelopes = load_envelopes(element, [model1, model2, model3])
```

## Shell Internal Forces

Stress resultants in shell elements.

```julia
forces = ShellInternalForces(shell, model)

# Membrane forces [N/m]
forces.Nxx, forces.Nyy, forces.Nxy

# Bending moments [N*m/m]  
forces.Mxx, forces.Myy, forces.Mxy

# Principal values
M1, M2, θ = principal_moments(forces)
σ_vm = von_mises_stress(forces, shell.thickness)
```

## Displacements

```julia
# Frame elements
disp = ElementDisplacements(element, model)
δ = disp.v(0.5)  # Transverse displacement at midspan

# Shell elements
disp = ShellDisplacements(shell, model)
w = max_deflection(disp)
```

---

# Tributary Area Computation

Asap includes algorithms for computing tributary areas—essential for distributing floor loads to beams.

## Straight Skeleton (Isotropic)

For two-way slab behavior:

```julia
using Meshes

vertices = [Point(0,0), Point(4,0), Point(4,3), Point(0,3)]
tributaries = get_tributary_polygons(vertices)

# Each edge gets a TributaryPolygon with parametric coordinates
for trib in tributaries
    area = trib.area
    coords = vertices(trib)  # Absolute coordinates
end
```

## One-Way Distribution

For one-way slabs with specified span direction:

```julia
tributaries = get_tributary_polygons(vertices; axis=[1.0, 0.0])  # Span in X
```

## Voronoi Tributaries

For column tributary areas:

```julia
column_tribs = compute_voronoi_tributaries(floor_polygon, column_points)
```

## Span Calculations

```julia
span_info = get_polygon_span(vertices)
Ln = short_span(span_info)
Ll = long_span(span_info)
```

---

# Force Density Method

Form-finding for cable and membrane structures.

```julia
# Define network
nodes = [
    FDMnode(0, 0, 0, false),  # Fixed
    FDMnode(10, 0, 0, false), # Fixed
    FDMnode(5, 5, 0, true),   # Free (will find form)
]

elements = [
    FDMelement(nodes, 1, 3, 1.0),  # q = force density
    FDMelement(nodes, 2, 3, 1.0),
]

loads = [FDMload(nodes[3], [0, 0, -100])]

network = Network(nodes, elements, loads)
solve!(network)

# Free node position after form-finding
nodes[3].position
```

---

# Acknowledgments

**Original Author**: [Keith Janghyun Lee](https://github.com/keithjlee) created the original Asap.jl package. This fork extends his excellent work with additional capabilities.

**Shell Formulations**: Based on [FinEtoolsFlexStructures.jl](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl) by Petr Krysl (MIT License). The T3FF shell element with DSG shear technology is adapted from his work.

**Dependencies**: This package builds on [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), and [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl).

---

# Citing

When using or extending this software for research purposes, please cite:

### Bibtex

```
@software{lee_2024_10581560,
  author       = {Lee, Keith Janghyun},
  title        = {Asap.jl},
  month        = jan,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v0.1},
  doi          = {10.5281/zenodo.10581560},
  url          = {https://doi.org/10.5281/zenodo.10581560}
}
```

### APA

Lee, K. J. (2024). Asap.jl (v0.1). Zenodo. https://doi.org/10.5281/zenodo.10581560

---

# Extensions

See also:

- [AsapToolkit](https://github.com/keithjlee/AsapToolkit) - Additional utilities and post-processing
- [AsapOptim](https://github.com/keithjlee/AsapOptim) - Structural optimization
- [AsapHarmonics](https://github.com/keithjlee/AsapHarmonics) - Harmonic analysis
