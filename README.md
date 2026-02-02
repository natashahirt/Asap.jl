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
- **Modal analysis** - Natural frequencies and mode shapes
- **Nonlinear analysis** - P-delta, linear buckling (frames + shells), Newton-Raphson pushover
- **Unified solve! API** - Single entry point with symbol dispatch: `solve!(model, :buckling)`
- **Mixed models** - Combined frame + shell structures in one model
- **Tributary areas** - Straight skeleton and Voronoi algorithms for load distribution
- **Area loads** - Unified surface pressure API for shells
- **Shell queries** - Spatial queries and region integration for design strips
- **Unit system** - Canonical source for structural engineering units (see below)

### Unit System (Canonical Source)

Asap is the **canonical source** for structural engineering units in the ecosystem. Other packages (`StructuralSizer`, `StructuralSynthesizer`) import and re-export units from here.

**US Customary Units:**
| Unit | Symbol | Definition |
|------|--------|------------|
| `kip` | kip | 1000 lbf (kilopound-force) |
| `ksi` | ksi | 1000 psi (kips per square inch) |
| `psf` | psf | lbf/ft² (pounds per square foot) |
| `ksf` | ksf | 1000 psf (kips per square foot) |
| `pcf` | pcf | lb/ft³ (pounds per cubic foot) |

**Type Aliases (for function signatures):**
- `Length`, `Area`, `Volume` - Geometric quantities
- `Pressure`, `Force`, `Moment`, `Torque` - Mechanical quantities
- `LinearLoad`, `Density`, `Acceleration` - Derived quantities
- `SecondMomentOfArea`, `TorsionalConstant`, `WarpingConstant` - Section properties

**Conversion Helpers:**
```julia
to_ksi(50u"MPa")      # → 7.252 (Float64 in ksi)
to_inches(1u"m")      # → 39.37 (Float64 in inches)
to_meters(12u"ft")    # → 3.658 (Float64 in meters)
to_kip(100u"kN")      # → 22.48 (Float64 in kip)
```

**Usage:**
```julia
using Asap
using Unitful

# Register units for u"..." string macro
Unitful.register(Asap)

# Now you can use:
load = 100u"psf"
stress = 50u"ksi"
length = 30u"ft"
```

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

# Modal analysis (unified API)
modal_result = solve!(model, :modal; n=10)
print_modal_summary(modal_result)

# Stability check (buckling + P-delta)
stab = solve!(model, :stability)
println("\\nBuckling load factor: ", round(stab.buckling.load_factors[1], digits=2))
println("P-delta amplification: ", round(stab.pdelta.amplification, digits=3))

# Access results
println("\\nFirst 3 natural frequencies: ", round.(modal_result.frequencies[1:3], digits=2), " Hz")
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
# Create materials inline
concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
steel = ShellMaterial(E=200u"GPa", ν=0.3, ρ=7850u"kg/m^3", name=:steel)
aluminum = ShellMaterial(E=70u"GPa", ν=0.33, ρ=2700u"kg/m^3", name=:aluminum)

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

# Or create sections directly
section_150mm = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
section_200mm = ShellSection(0.20u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
```

### `Shell()` - Create Shell Elements

Automatically triangulate any polygon into `ShellTri3` elements. Works with triangles, quads, or any convex/concave polygon.

```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")

# Default refinement (n=4)
shells = Shell((n1, n2, n3, n4), section)

# Finer mesh for higher accuracy
shells = Shell((n1, n2, n3, n4), 8, section)

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

Asap provides a **unified `solve!` API** with symbol dispatch for all analysis types.

## Unified API

```julia
# Default: linear static
solve!(model)                          # K·u = F

# Explicit analysis types
solve!(model, :static)                 # Linear static
solve!(model, :modal; n=6)             # Modal analysis  
solve!(model, :buckling; n=4)          # Linear buckling
solve!(model, :pdelta; max_iter=10)    # P-delta iteration
solve!(model, :nonlinear; n_steps=10)  # Newton-Raphson pushover

# Composite analyses (run multiple at once)
solve!(model, :stability)              # Buckling + P-delta → NamedTuple
solve!(model, :all)                    # All valid for this model → NamedTuple
```

**Model-aware dispatch:** `solve!(model, :all)` only runs analyses valid for the model type.

```julia
available_analyses(frame_model)  # → (:static, :modal, :buckling, :pdelta, :nonlinear)
available_analyses(shell_model)  # → (:static, :modal, :buckling)
available_analyses(truss_model)  # → (:static, :buckling)
```

## Static Analysis

```julia
solve!(model)  # or solve!(model, :static)

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
result = solve!(model, :modal)       # 6 modes (default)
result = solve!(model, :modal; n=10) # 10 modes

# Access results
result.frequencies   # Natural frequencies [Hz]
result.periods       # Periods [s]
result.mode_shapes   # Mode shape matrix (columns = modes)
result.omegas        # Angular frequencies [rad/s]

# Pretty print summary
print_modal_summary(result)

# Mass matrix options
result = solve!(model, :modal; n=6, mass_type=MASS_LUMPED)

# Legacy API still works
result = modal!(model; n=10)
```

Mass types: `MASS_CONSISTENT` (default), `MASS_LUMPED`, `MASS_CONSISTENT_NO_ROTATION`, `MASS_LUMPED_NO_ROTATION`

## Buckling Analysis

Linear buckling eigenvalue problem: (K + λ·Kg)φ = 0

```julia
result = solve!(model, :buckling; n=4)

# Access results
result.load_factors  # Critical load multipliers [λ₁, λ₂, ...]
result.mode_shapes   # Buckling mode shapes
result.n_modes       # Number of modes computed

# Interpretation
# λ > 1.0 → Structure is stable under current loads
# λ < 1.0 → Structure will buckle before reaching full load
# λ = 1.0 → At critical load

print_buckling_summary(result)

# Convenience functions
critical_load_factor(model)  # Just λ₁
is_stable(model)             # λ₁ > 1.0?
```

**Example: Column buckling check**

```julia
# Cantilever column under axial load
nodes = [
    Node([0u"m", 0u"m", 0u"m"], :fixed),
    Node([0u"m", 0u"m", 5u"m"], :free)
]
column = Element(nodes[1], nodes[2], section)
load = NodeForce(nodes[2], [0u"N", 0u"N", -500u"kN"])

model = Model(nodes, [column], [load])
result = solve!(model, :buckling)

if result.load_factors[1] > 1.5
    println("Column OK with safety factor $(result.load_factors[1])")
else
    println("WARNING: Near buckling!")
end
```

## P-Delta Analysis

Second-order effects from gravity loads acting through lateral displacements.

```julia
result = solve!(model, :pdelta; max_iter=10, tol=1e-3)

# Access results
result.converged      # Did iteration converge?
result.iterations     # Number of iterations
result.amplification  # Second-order amplification factor (B₂)
result.max_drift_ratio

print_pdelta_summary(result)

# AISC B₂ factor directly
B2 = B2_factor(model)
```

**When to use P-delta:**
- Tall buildings (>10 stories)
- Slender frames
- When B₂ > 1.1 per AISC 360

**Example: Frame stability check**

```julia
model = Model(nodes, columns ∪ beams, gravity_loads ∪ lateral_loads)
result = solve!(model, :pdelta)

if result.amplification > 1.5
    println("WARNING: Significant P-delta effects")
    println("Amplification factor: $(round(result.amplification, digits=2))")
end
```

## Nonlinear Static Analysis

Full Newton-Raphson with incremental loading—for pushover analysis and capacity curves.

```julia
result = solve!(model, :nonlinear; n_steps=20, max_iter=20, tol=1e-4)

# Access results
result.converged           # Did all steps converge?
result.load_factors        # λ at each step [0, 0.05, 0.10, ...]
result.displacements       # u at each step
result.iterations_per_step # Iterations needed
result.equilibrium_error   # Final residual at each step

print_nonlinear_summary(result)

# Extract capacity curve for plotting
curve = capacity_curve(result)
# curve.displacement  # Max displacement at each step
# curve.load_factor   # Load level at each step

# Alias for seismic assessment
result = pushover!(model; n_steps=20)
```

**Example: Pushover capacity curve**

```julia
# Apply lateral load pattern
model = Model(nodes, elements, [lateral_load_at_roof])
result = solve!(model, :nonlinear; n_steps=20)

curve = capacity_curve(result)
# Plot: curve.displacement vs curve.load_factor
```

## Stability Check (Composite)

Run buckling + P-delta together for comprehensive stability assessment.

```julia
results = solve!(model, :stability)

# Access both results
results.buckling.load_factors[1]  # Critical buckling λ
results.pdelta.amplification      # P-delta amplification

# Quick summary
println("Buckling λ₁ = ", results.buckling.load_factors[1])
println("P-delta B₂ = ", results.pdelta.amplification)
```

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

# Principal values (returns NamedTuple)
pm = principal_moments(forces)  # (M1=..., M2=..., θ=...)
pm.M1, pm.M2, pm.θ              # Named access
M1, M2, θ = principal_moments(forces)  # Destructuring also works

pf = principal_forces(forces)   # (N1=..., N2=..., θ=...)
σ_vm = von_mises_stress(forces, shell.thickness)
```

## Shell Queries & Region Integration

Query shell results by location or integrate over regions—essential for design strip calculations.

### Spatial Queries

```julia
# Find shell triangles at a point (returns all touching triangles)
tris = shell_tris_at_point(model, (5.0, 3.0); tol=0.01)

# Find triangles whose centroids fall in a polygon
column_strip = [(0.0, 0.0), (10.0, 0.0), (10.0, 2.5), (0.0, 2.5)]
tris = shell_tris_in_region(model, column_strip)

# Get element centroid (returns NamedTuple)
c = shell_centroid(shell)      # (x=..., y=...)
c.x, c.y                       # Named access
cx, cy = shell_centroid(shell) # Destructuring also works

c3d = shell_centroid_3d(shell) # (x=..., y=..., z=...)
```

### Bending Moments

Multiple dispatch patterns for querying moments:

```julia
# Single element
M = bending_moments(shell, model)  # Returns [Mxx, Myy, Mxy]

# Multiple elements
Ms = bending_moments(shells, model)  # Returns Vector{[Mxx, Myy, Mxy]}

# At a point (averages if multiple triangles touch the point)
M = bending_moments(model, (5.0, 3.0); tol=0.01)

# Integrate over a region (polygon)
strip = [(0.0, 0.0), (10.0, 0.0), (10.0, 2.5), (0.0, 2.5)]
result = bending_moments(model; polygon=strip)
# Returns NamedTuple:
#   result.Mxx, .Myy, .Mxy         - Area-weighted totals
#   result.Mxx_avg, .Myy_avg, ...  - Average moment intensity (N·m/m)
#   result.Mxx_max, .Myy_max, ...  - Peak absolute values
#   result.area                    - Total integrated area (m²)
#   result.shell_tris              - Shells in the region

# Query multiple points
pts = [(1.0, 1.0), (5.0, 5.0), (9.0, 9.0)]
moments = bending_moments(model; pts=pts, tol=0.01)

# From a filtered element set (e.g., specific slab)
slab_tris = shell_tris_in_region(model, slab_boundary)
result = bending_moments(slab_tris, model; polygon=strip_boundary)
```

### Design Strip Workflow

Typical workflow for ACI-style design strip analysis:

```julia
# 1. Define strip geometry (column strip, middle strip, etc.)
column_strip = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]

# 2. Integrate FEA moments over strip region
result = bending_moments(model; polygon=column_strip)

# 3. Use results for reinforcement design
Mu_design = result.Mxx_max  # Peak moment for strength design
Mu_avg = result.Mxx_avg     # Average for minimum steel calculation
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

# Acknowledgments

**Original Author**: [Keith Janghyun Lee](https://github.com/keithjlee) created the original Asap.jl package. This fork extends his excellent work with additional capabilities.

**Shell Formulations**: Based on [FinEtoolsFlexStructures.jl](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl) by Petr Krysl (MIT License). The T3FF shell element with DSG shear technology is adapted from his work.

**Shell Buckling**: The shell geometric stiffness implementation follows classical plate theory (Cook et al., Przemieniecki). The [FEniCSx-Shells](https://github.com/FEniCS-Shells/fenicsx-shells) nonlinear Naghdi demo was a helpful reference for understanding shell strain measures.

**Dependencies**: This package builds on [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), and [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl).

---

# Extensions

See also:

- [AsapToolkit](https://github.com/keithjlee/AsapToolkit) - Additional utilities and post-processing
- [AsapOptim](https://github.com/keithjlee/AsapOptim) - Structural optimization
- [AsapHarmonics](https://github.com/keithjlee/AsapHarmonics) - Harmonic analysis
