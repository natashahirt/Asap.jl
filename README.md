[![DOI](https://zenodo.org/badge/426740094.svg)](https://zenodo.org/doi/10.5281/zenodo.10581559)

![](READMEassets/forces-axo.png)

# Asap.jl

Asap is...

- Another Structural Analysis Package
- results As Soon As Possible
- Analysis of Structures Avec Programming
- A Simple Analysis Please

Designed first-and-foremost for information-rich data structures and ease of querying, but always with performance in mind.

See also: [AsapToolkit](https://github.com/keithjlee/AsapToolkit), [AsapOptim](https://github.com/keithjlee/AsapOptim), [AsapHarmonics](https://github.com/keithjlee/AsapHarmonics).

# Installation

Asap.jl is now a registered Julia package. Install through package mode in the REPL:

```julia
pkg> add Asap
```

or

```julia
using Pkg
Pkg.Add("Asap")
```

## Citing

When using or extending this software for research purposes, please cite using the following:

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

### Other styles

Or find a pre-written citation in the style of your choice [here](https://zenodo.org/records/10724610) (see the Citation box on the right side). E.g., for APA:

```
Lee, K. J. (2024). Asap.jl (v0.1). Zenodo. https://doi.org/10.5281/zenodo.10581560
```

# Extensions, Related packages

See [AsapToolkit.jl](https://github.com/keithjlee/AsapToolkit) for even more utility and post-processing functions.

# Units (Unitful.jl)

Asap uses [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) for type-safe unit handling. All physical quantities (distances, forces, pressures, etc.) should be specified with units:

```julia
using Unitful

# Node positions with units
n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
n2 = Node([10.0u"ft", 0.0u"ft", 0.0u"ft"], :free)  # Can mix unit systems

# Section properties
sec = Section(
    11.8u"inch^2",      # Area
    29e3u"ksi",         # Young's modulus
    11.5e3u"ksi",       # Shear modulus
    310.0u"inch^4",     # Ix
    100.0u"inch^4",     # Iy
    15.0u"inch^4"       # J
)

# Loads with units
load = NodeForce(n2, [0.0u"kN", -50.0u"kN", 0.0u"kN"])
line_load = LineLoad(element, [0.0u"kip/ft", -1.5u"kip/ft", 0.0u"kip/ft"])
```

Asap automatically converts all quantities to SI units internally (meters, Newtons, Pascals) and stores results in SI. When retrieving results, you can convert back to your preferred units:

```julia
# Displacement in inches
d_inches = uconvert(u"inch", node.displacement[1])

# Reaction in kips
R_kips = uconvert(u"kip", node.reaction[2])
```

# Usage

A structural model is defined by:

```julia
model = Model(nodes, elements, loads)
```

and solved via:

```julia
solve!(model)
```

Which finds the unknown nodal displacement field, $u = S^{-1}(P-P_f)$ where:

- $S$ is the global stiffness matrix (often called $K$)
- $P$ is the global external load vector
- $P_f$ is the fixed end forces induced by element loads (such as a line load on an element)

## `Node`

```julia
mutable struct Node <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Symbol
end
```

We begin with the primary information carrier for structural analysis: nodes with *n* independent degrees of freedom (DOF). They are defined by a spatial position in $\mathbb{R}^3$ as well as a vector of booleans that indicate which DOFs are free to move under load, in order: $T_x, T_y, T_z, R_x, R_y, R_z$ where $T$ is a translational DOF and $R$ is a rotational DOF. E.g.:

```julia
node = Node([0u"m", 15.5u"m", 12.0u"m"], [false, false, false, true, true, true])
```

This defines a node at $x = 0; y = 15.5; z = 12$ meters with a *pinned* support (i.e., translational DOFs are fixed, but rotational DOFs are not).

Some common boundary conditions are provided to you as symbols to use in the constructor:

```julia
node = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)    # all free DOFs
node = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)   # all fixed DOFs
node = Node([0.0u"m", 0.0u"m", 0.0u"m"], :xfree)   # all DOFs are fixed except Tx
node = Node([0.0u"m", 0.0u"m", 0.0u"m"], :xfixed)  # all DOFs are free except Tx
```

Nodes can also include an optional identifier represented as a symbol:

```julia
pin_support = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned, :pinsupport)
roller_support = Node([15.1u"m", 0.0u"m", 0.0u"m"], :xfree, :rollersupport)
free_nodes = [Node([rand()*10u"m", rand()*10u"m", rand()*10u"m"], :free, :freenodes) for _ = 1:10]

nodes = [pin_support; roller_support; free_nodes]
```

This allows you to index into a vector of nodes using the identifier:

```julia
all_free_nodes = nodes[:freenodes] #returns a vector of nodes with the :freenode identifier
```

Or find the indices of nodes in a vector of nodes that have a given id:

```
i_free_nodes = findall(nodes, :freenodes)
```

## `Element`

```julia
mutable struct Element{R<:Release} <: FrameElement{R}
    section::Section #cross section
    nodeStart::Node #start node
    nodeEnd::Node #end position
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    Q::Vector{Float64} # fixed end forces in GCS
    R::Matrix{Float64} # transformation matrix
    Ψ::Float64 #roll angle
    LCS::Vector{Vector{Float64}} #local coordinate frame (X, y, z)
    forces::Vector{Float64} #elemental forces in LCS
    id::Symbol #optional identifier
end
```

Elements are defined by their start and end nodes, a cross-section, and an optional identifier.

### `Section`

A section defines the mechanical and material properties of an element:

```julia
ibeam_section = Section(A, E, G, Ix, Iy, J)
```

where:

- `A`: area (e.g., `u"m^2"`, `u"inch^2"`)
- `E`: Young's Modulus (e.g., `u"Pa"`, `u"GPa"`, `u"ksi"`)
- `G`: Shear Modulus (e.g., `u"Pa"`, `u"GPa"`, `u"ksi"`)
- `Ix`: Moment of inertia in strong axis (e.g., `u"m^4"`, `u"inch^4"`)
- `Iy`: Moment of inertia in weak axis
- `J`: Torsional constant

All values must have Unitful units. Asap automatically converts to SI internally.

Example:

```julia
# SI units
steel_section = Section(
    0.01u"m^2",      # A
    200.0u"GPa",     # E
    77.0u"GPa",      # G
    1e-4u"m^4",      # Ix
    5e-5u"m^4",      # Iy
    1e-5u"m^4"       # J
)

# Imperial units
w_section = Section(
    11.8u"inch^2",   # A
    29e3u"ksi",      # E
    11.5e3u"ksi",    # G
    310.0u"inch^4",  # Ix
    37.1u"inch^4",   # Iy
    1.39u"inch^4"    # J
)
```

You can then define an element via:

```julia
element = Element(pin_support, roller_support, steel_section)
element_with_id = Element(pin_support, rand(free_nodes), steel_section, :randomelement)
```

### `Element` roll axis

Elements have a default roll angle with respect to its longitudinal axis of $\pi/2$, which corresponds to keeping the strong bending axis flat against the XY plane. If you wish to change this, you can change it by accessing the Ψ parameter:

```julia
element.Ψ = pi
```

### `Element` release

Elements can have partial DOF releases to decouple nodal displacements from element end displacements. This is the process of adding hinges to one or both ends of the beam. You can do this in the construction of an element through the optional argument `release`:

```julia
released_element = Element(pin_support, roller_support, ibeam_section; release = :freefixed)
```

By default, no releases are performed (i.e., `release = :fixedfixed`). You can choose between:

- `:freefixed` create a hinge in the beginning node of the element
- `:fixedfree` create a hinge in the ending node of the element
- `:freefree` create hinges on both ends of the element
- `:joist` create hinges on both ends of the element with the exception of torsional DOFs.

# Loads

Loads can be applied to nodes and elements. All load values require Unitful units.

## `NodeForce`

A `NodeForce` is defined on a node with a force vector [Fx, Fy, Fz]:

```julia
load1 = NodeForce(free_nodes[1], [0.0u"kN", 0.0u"kN", -150.0u"kN"])

# Or in imperial units
load2 = NodeForce(free_nodes[2], [0.0u"lbf", 0.0u"lbf", -10.0u"kip"])
```

## `NodeMoment`

A `NodeMoment` is defined on a node with a moment vector [Mx, My, Mz]:

```julia
moment_load = NodeMoment(free_nodes[5], [40.0u"kN*m", 0.0u"kN*m", 0.0u"kN*m"])

# Or in imperial
moment_load = NodeMoment(free_nodes[5], [100.0u"kip*ft", 0.0u"kip*ft", 0.0u"kip*ft"])
```

## `LineLoad`

A `LineLoad` is a distributed load [wx, wy, wz] in force/length applied along an element:

```julia
# 10 kN/m downward load
snow_load = LineLoad(element, [0.0u"kN/m", 0.0u"kN/m", -10.0u"kN/m"])

# 1.5 kip/ft in imperial
floor_load = LineLoad(element, [0.0u"kip/ft", 0.0u"kip/ft", -1.5u"kip/ft"])
```

## `PointLoad`

A `PointLoad` is applied at a normalized position $0<x<1$ along an element:

```julia
# 20 kN lateral load at quarter point
sideways_load = PointLoad(element, 0.25, [20.0u"kN", 0.0u"kN", 0.0u"kN"])

# 50 kip vertical at midspan
midspan_load = PointLoad(element, 0.5, [0.0u"kip", 0.0u"kip", -50.0u"kip"])
```

## `TributaryLoad`

A `TributaryLoad` represents a piecewise-linear distributed load, typically derived from a tributary polygon. This is useful for floor loads where the tributary area varies along the beam length.

```julia
mutable struct TributaryLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    positions::Vector{Float64}       # Breakpoints [0, s1, s2, ..., 1], normalized
    widths::Vector{Length}           # Tributary widths at each position (Unitful)
    pressure::QuantityPressure       # Load intensity (Pa)
    direction::NTuple{3, Float64}    # Unit load direction
    loadID::Int64
    id::Symbol
end
```

The load is defined by:

- `positions`: Normalized positions along the beam (0 = start, 1 = end). Must be sorted and within [0, 1].
- `widths`: Tributary widths at each position with distance units (e.g., `u"m"`, `u"ft"`, `u"inch"`). Converted to meters internally.
- `pressure`: Surface pressure (force/area), e.g., `5000.0u"Pa"` or `100.0u"psf"` for floor loads.
- `direction`: Unit vector for load direction. Default is `(0.0, 0.0, -1.0)` (gravity/downward).

### Example: Uniform tributary load

A beam with uniform 2m tributary width and 5 kPa floor load:

```julia
trib_load = TributaryLoad(
    element,
    [0.0, 1.0],                    # Full span
    [2.0u"m", 2.0u"m"],            # 2m width everywhere
    5000.0u"Pa",                   # 5 kPa
    (0.0, 0.0, -1.0)               # Downward
)
```

This is equivalent to a `LineLoad` of $w = 2.0 \times 5000 = 10{,}000$ N/m.

### Example: Imperial units

The same load in imperial units:

```julia
trib_load = TributaryLoad(
    element,
    [0.0, 1.0],
    [6.56u"ft", 6.56u"ft"],        # ~2m in feet
    104.4u"psf"                     # ~5 kPa in psf
)
```

### Example: Varying tributary width

A beam where tributary width varies from 1m at start to 3m at end:

```julia
trib_load = TributaryLoad(
    element,
    [0.0, 1.0],                    # Full span
    [1.0u"m", 3.0u"m"],            # 1m → 3m
    5000.0u"Pa"
)
```

### Example: Piecewise tributary from polygon

For complex tributary shapes (e.g., from a straight skeleton algorithm):

```julia
# Tributary polygon with multiple breakpoints
positions = [0.0, 0.2, 0.5, 0.8, 1.0]
widths = [0.5u"m", 1.5u"m", 2.0u"m", 1.5u"m", 0.5u"m"]  # Peak width at midspan

trib_load = TributaryLoad(element, positions, widths, 4000.0u"Pa")
```

### Computing line load intensities

The `intensities` function returns the equivalent line load (N/m) at each breakpoint:

```julia
w = intensities(trib_load)  # Returns widths .* pressure in N/m
```

This is useful for debugging or comparing with `LineLoad` values.

### Updating pressure

The `pressure` field is mutable, allowing you to update the load magnitude without recreating the load object:

```julia
# Change floor type (e.g., from office to storage)
trib_load.pressure = 7500.0u"Pa"
solve!(model; reprocess = true)
```

## `Model`

```julia
mutable struct Model{E,L} <: AbstractModel
    nodes::Vector{Node}
    elements::Vector{E}
    loads::Vector{L}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool} #vector of DOFs
    nDOFs::Int64
    freeDOFs::Vector{Int64} #free DOF indices
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64} # global stiffness
    P::Vector{Float64} # external loads
    Pf::Vector{Float64} # element end forces
    u::Vector{Float64} # nodal displacements
    reactions::Vector{Float64} # reaction forces
    compliance::Float64 #structural compliance
    tol::Float64
    processed::Bool
end
```

A model is assembled from a collection of nodes, elements, and loads:

```julia
nodes = [pin_support; roller_support; free_nodes]
elements = [element, released_element]
loads = [load1, snow_load, sideways_load]

model = Model(nodes, elements, loads)
```

## Solving

The primary unknown field we are trying to find is `u`, the vector of all nodal DOFs (in order of assembly) in which equilibrium holds. We can find this via:

```julia
solve!(model)
```

You can access the solved field via:

```julia
u = model.u
```

Or directly from the populated fields in the nodes (results are Unitful quantities):

```julia
# Displacement is a vector of Unitful quantities [Tx, Ty, Tz, Rx, Ry, Rz]
node2_displacement = model.nodes[2].displacement
dy = node2_displacement[2]  # Returns e.g. -0.0023 m

# Convert to your preferred units
dy_mm = uconvert(u"mm", dy)   # -2.3 mm
dy_inch = uconvert(u"inch", dy)  # -0.0906 inch
```

If a node has a restrained DOF, you can find its reaction from:

```julia
# Reaction is a vector of Unitful quantities [Fx, Fy, Fz, Mx, My, Mz]
roller_reaction = roller_support.reaction
Fy = roller_reaction[2]  # Returns e.g. 75000 N

# Convert to kN or kip
Fy_kN = uconvert(u"kN", Fy)    # 75 kN
Fy_kip = uconvert(u"kip", Fy)  # 16.86 kip
```

You can also find the end forces acting on an element via:

```julia
element_forces = model.elements[2].forces
```

Which gives a vector: $F_{x1}, F_{y1}, F_{z1}, M_{x1}, M_{y1}, M_{z1}, F_{x2}, F_{y2}, F_{z2}, M_{x2}, M_{y2}, M_{z2}$ where $1$ is the starting node and $2$ is the ending node, with all values defined in the *local coordinate system* of the beam.

## New loads

If you have a new set of loads, directly get the corresponding displacement via:

```julia
u_new = solve!(model, new_loads)
```

Or replace the vector of loads associated with the model and solve in place via:

```julia
solve!(model, new_loads)
```

## Updating values

If you change a value, such as the position of a node, reprocess the fields before solving by:

```julia
# Modify position (note: position is stored in meters)
model.nodes[2].position .+= [5.0u"m", 0.0u"m", 0.0u"m"]
solve!(model; reprocess = true)
```

# Trusses

For truss structures, with only 3 translational DOFs per node, there are separate data structures for `TrussNode`, `TrussElement`, `TrussSection`, and `TrussModel`, which can be defined similarly as above except:

1. `TrussNode`s are constructed with Unitful positions and a length-3 vector of booleans: `TrussNode([1.0u"m", 2.0u"m", 0.0u"m"], [true, true, false])`.
2. `TrussSection`s only require the area and modulus with units: `TrussSection(0.01u"m^2", 200u"GPa")`.
3. `TrussElement`s do not have releases or roll angles. By definition they are equivalent to `Element(...; release = :freefree)`
4. Only `NodeForce`s (with Unitful values) can be applied as loads for `TrussModel`s.

## `TrussNode`

```julia
mutable struct TrussNode <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Symbol
end
```

## `TrussElement`

```julia
mutable struct TrussElement <: AbstractElement
    section::Union{TrussSection,Section} #cross section
    nodeStart::TrussNode #start position
    nodeEnd::TrussNode #end position
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    R::Matrix{Float64} # transformation matrix
    forces::Vector{Float64} #elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing} #optional identifier
end
```

## `TrussModel`

```julia
mutable struct TrussModel <: AbstractModel
    nodes::Vector{TrussNode}
    elements::Vector{TrussElement}
    loads::Vector{NodeForce}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool} #vector of DOFs
    nDOFs::Int64
    freeDOFs::Vector{Int64} #free DOF indices
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64} # global stiffness
    P::Vector{Float64} # external loads
    u::Vector{Float64} # nodal displacements
    reactions::Vector{Float64} # reaction forces
    compliance::Float64 #structural compliance
    tol::Float64
    processed::Bool
end
```
