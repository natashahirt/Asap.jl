"""
Define loads applied to structure/elements
"""
abstract type AbstractLoad end
abstract type NodeLoad <: AbstractLoad end
abstract type ElementLoad{R} <: AbstractLoad end

"""
    NodeForce(node::AbstractNode, value::Vector{Quantity})

A force vector [Fx, Fy, Fz] in the global coordinate system applied to a node.
"""
mutable struct NodeForce <: NodeLoad
    node::Union{Node, TrussNode}
    value::Vector{QuantityForce}
    loadID::Int64
    id::Symbol
    
    function NodeForce(node::AbstractNode, value::Vector{<:Quantity}, id::Symbol = :force)
        @assert length(value) == 3 "load vector must be in R³ (GCS)"
        
        # Convert all to Newtons
        value_si = [uconvert(u"N", v) for v in value]
        force = new(node, value_si, 0, id)
        return force
    end

    function NodeForce(nodes::Vector{<:AbstractNode}, index::Integer, value::Vector{<:Quantity}, id::Symbol = :force)
        return NodeForce(nodes[index], value, id)
    end
end

"""
    NodeMoment(node::Node, value::Vector{Quantity})

A moment vector [Mx, My, Mz] in the global coordinate system applied to a node with rotational DOFs.
"""
mutable struct NodeMoment <: NodeLoad
    node::Node
    value::Vector{Moment}
    loadID::Int64
    id::Symbol
    
    function NodeMoment(node::Node, value::Vector{<:Quantity}, id::Symbol = :moment)
        @assert length(value) == 3 "Moment vector must be in R³ (GCS)"
        
        # Convert all to N*m
        value_si = [uconvert(u"N*m", v) for v in value]
        force = new(node, value_si, 0, id)
        return force
    end
end

"""
    LineLoad(element::Element, value::Vector{Quantity})

A distributed line load [wx, wy, wz] in (force/length) applied along an element in the global coordinate system.
"""
mutable struct LineLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    value::Vector{ForcePerLength}
    loadID::Int64
    id::Symbol

    function LineLoad(element::FrameElement{R}, value::Vector{<:Quantity}, id::Symbol = :lineload) where R
        @assert length(value) == 3 "load vector must be in R³ (GCS)"
        
        # Convert all to N/m
        value_si = [uconvert(u"N/m", v) for v in value]
        force = new{R}(element, value_si, 0, id)
        return force
    end
end

"""
    GravityLoad(element::Element, factor::Quantity)

A gravity load (negative global Z) applied along a member. 

Generates distributed load w = element.section.A * element.section.ρ * factor, where factor should be the appropriate acceleration due to gravity.
"""
mutable struct GravityLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    factor::QuantityAcceleration
    loadID::Int64
    id::Symbol

    function GravityLoad(element::FrameElement{R}, factor::Quantity, id::Symbol = :gravityload) where R
        # Convert to m/s²
        factor_si = uconvert(u"m/s^2", factor)
        force = new{R}(element, factor_si, 0, id)
        return force
    end
end


"""
    PointLoad(element::Element, position::Float64, value::Vector{Quantity})

A point load [Px, Py, Pz] applied in the global coordinate system at a distance `position` × `element.length` from the starting node.
`position` is dimensionless (fraction 0-1).
"""
mutable struct PointLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    position::Float64  # Dimensionless fraction
    value::Vector{QuantityForce}
    loadID::Int64
    id::Symbol

    function PointLoad(element::FrameElement{R}, position::Float64, value::Vector{<:Quantity}, id::Symbol = :pointload) where R
        @assert 0 < position < 1 "position must be ∈ ]0, 1["
        @assert length(value) == 3 "load vector must be in R³ (GCS)"
        
        # Convert all to Newtons
        value_si = [uconvert(u"N", v) for v in value]
        force = new{R}(element, position, value_si, 0, id)
        return force
    end
end

"""
    TributaryLoad(element, positions, widths, pressure, direction)

A piecewise-linear distributed load from a tributary polygon.

The load intensity at each position is `widths[i] * pressure`, varying linearly between 
breakpoints. This enables exact fixed-end force computation without discretization.

## Arguments
- `element::FrameElement`: The beam element
- `positions::Vector{Float64}`: Breakpoint positions along beam, normalized [0, 1]
- `widths::Vector{<:Quantity}`: Tributary widths at each position (distance units, e.g. m, ft)
- `pressure::Quantity{Pressure}`: Load intensity (Pa or N/m²)
- `direction::NTuple{3, Float64}`: Unit load direction vector (default gravity: (0,0,-1))

## Updating
The `pressure` field is mutable - update it and re-solve to change load magnitude
while keeping the same tributary geometry.
"""
mutable struct TributaryLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    positions::Vector{Float64}       # Breakpoints [0, s1, s2, ..., 1], normalized
    widths::Vector{Length}           # Tributary widths at each position (m)
    pressure::QuantityPressure       # Load intensity (Pa)
    direction::NTuple{3, Float64}    # Unit load direction
    loadID::Int64
    id::Symbol

    function TributaryLoad(
        element::FrameElement{R},
        positions::Vector{Float64},
        widths::Vector{<:Quantity},
        pressure::Quantity,
        direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
        id::Symbol = :tributaryload
    ) where R
        n = length(positions)
        @assert n >= 2 "positions must have at least 2 elements"
        @assert length(widths) == n "widths must have same length as positions"
        @assert all(0.0 .<= positions .<= 1.0) "positions must be in [0, 1]"
        @assert issorted(positions) "positions must be sorted"
        
        # Convert widths to meters
        widths_si = [uconvert(u"m", w) for w in widths]
        @assert all(ustrip.(widths_si) .>= 0.0) "widths must be non-negative"
        
        # Normalize direction
        dir_len = sqrt(sum(direction .^ 2))
        dir_norm = dir_len > 1e-12 ? direction ./ dir_len : (0.0, 0.0, -1.0)
        
        # Convert pressure to Pa
        pressure_si = uconvert(u"Pa", pressure)
        
        load = new{R}(element, positions, widths_si, pressure_si, dir_norm, 0, id)
        return load
    end
end

"""Compute line load intensities (N/m) at each breakpoint."""
function intensities(load::TributaryLoad)::Vector{Float64}
    p = to_pascals(load.pressure)
    w = [to_meters(width) for width in load.widths]
    return w .* p  # width (m) × pressure (N/m²) = N/m
end

# =============================================================================
# Area Loads (Pressure on Shell Elements)
# =============================================================================

"""
    AreaLoad(shells, pressure; distribute_to=:nodes, interior_beams=[], axis=nothing, direction=(0,0,-1))

A uniform pressure load over shell element surfaces.

Completes the load type progression: `PointLoad` → `LineLoad` → `AreaLoad`

# Arguments
- `shells`: Shell elements defining the loaded area
  - Single element: `ShellTri3`
  - Multiple elements: `Vector{<:ShellElement}`
- `pressure::Quantity{Pressure}`: Load intensity (Pa or N/m²)

# Keywords
- `distribute_to`: How load is distributed
  - `:nodes` (default): FEM approach - equivalent nodal forces on shell nodes
  - `Vector{FrameElement}`: Tributary approach - distribute to supporting beams/columns
- `interior_beams`: Interior beams that split the slab region (for tributary distribution)
  - These beams receive tributary load from both sides
  - Only used when `distribute_to` is a beam vector
- `axis`: For tributary distribution, controls one-way vs two-way
  - `nothing` (default): Two-way isotropic (straight skeleton)
  - `(1.0, 0.0)`: One-way along X (load to Y-perpendicular edges)
  - `(0.0, 1.0)`: One-way along Y (load to X-perpendicular edges)
- `direction::NTuple{3, Float64}`: Load direction vector (default: gravity)

# Examples
```julia
# FEM approach: forces at shell nodes (for shell analysis)
load = AreaLoad(shells, 5000u"Pa")

# Tributary to edge beams only
load = AreaLoad(shells, 5000u"Pa"; distribute_to=edge_beams)

# Tributary to edge + interior beams (interior beams get load from both sides)
load = AreaLoad(shells, 5000u"Pa"; 
    distribute_to=vcat(edge_beams, interior_beam),
    interior_beams=[interior_beam]
)

# One-way slab spanning along X
load = AreaLoad(shells, 5000u"Pa"; distribute_to=beams, axis=(1.0, 0.0))
```
"""
mutable struct AreaLoad <: AbstractLoad
    shells::Vector{<:ShellElement}
    pressure::QuantityPressure
    direction::NTuple{3, Float64}
    distribute_to::Union{Symbol, Vector{<:FrameElement}}
    interior_beams::Vector{<:FrameElement}
    axis::Union{Nothing, NTuple{2, Float64}}
    loadID::Int64
    id::Symbol
    
    # Internal: cached tributary loads (computed lazily during process!)
    _tributary_loads::Union{Nothing, Vector{TributaryLoad}}

    function AreaLoad(
        shells::Vector{S},
        pressure::Quantity;
        distribute_to::Union{Symbol, Vector{<:FrameElement}} = :nodes,
        interior_beams::Vector{<:FrameElement} = FrameElement[],
        axis::Union{Nothing, NTuple{2, <:Real}, Vector{<:Real}} = nothing,
        direction::NTuple{3, Float64} = (0.0, 0.0, -1.0),
        id::Symbol = :areaload
    ) where S <: ShellElement
        # Normalize direction
        dir_len = sqrt(sum(direction .^ 2))
        dir_norm = dir_len > 1e-12 ? direction ./ dir_len : (0.0, 0.0, -1.0)
        
        # Convert pressure to Pa
        pressure_si = uconvert(u"Pa", pressure)
        
        # Normalize axis to tuple
        axis_tuple = if isnothing(axis)
            nothing
        elseif axis isa Vector
            (Float64(axis[1]), Float64(axis[2]))
        else
            (Float64(axis[1]), Float64(axis[2]))
        end
        
        new(shells, pressure_si, dir_norm, distribute_to, interior_beams, axis_tuple, 0, id, nothing)
    end
end

# Convenience: single shell element
function AreaLoad(
    shell::ShellElement,
    pressure::Quantity;
    kwargs...
)
    AreaLoad([shell], pressure; kwargs...)
end

"""
    AreaLoad(shells, load_vector; ...)

Create an area load from a 3D load vector (pressure × direction combined).

This is a convenience interface where the load vector encodes both magnitude
and direction: `[px, py, pz]` where the magnitude is `sqrt(px² + py² + pz²)`
and direction is the normalized vector.

# Examples
```julia
# 90 Pa downward (Z-)
load = AreaLoad(shells, [0, 0, -90]u"Pa")

# 50 Pa at 45° in XZ plane
load = AreaLoad(shells, [35.36, 0, -35.36]u"Pa")

# Using N/m² (same as Pa)
load = AreaLoad(shells, [0, 0, -5000]u"N/m^2")
```

This is equivalent to:
```julia
magnitude = norm([px, py, pz])
direction = [px, py, pz] / magnitude
AreaLoad(shells, magnitude; direction=Tuple(direction))
```
"""
function AreaLoad(
    shells::Vector{S},
    load_vector::AbstractVector{<:Quantity};
    distribute_to::Union{Symbol, Vector{<:FrameElement}} = :nodes,
    interior_beams::Vector{<:FrameElement} = FrameElement[],
    axis::Union{Nothing, NTuple{2, <:Real}, Vector{<:Real}} = nothing,
    id::Symbol = :areaload
) where S <: ShellElement
    length(load_vector) == 3 || error("Load vector must have 3 components [px, py, pz]")
    
    # Convert to Pa and extract components
    px = ustrip(u"Pa", load_vector[1])
    py = ustrip(u"Pa", load_vector[2])
    pz = ustrip(u"Pa", load_vector[3])
    
    # Compute magnitude and direction
    magnitude = sqrt(px^2 + py^2 + pz^2)
    
    if magnitude < 1e-12
        # Zero load - use default direction
        direction = (0.0, 0.0, -1.0)
        pressure = 0.0u"Pa"
    else
        direction = (px / magnitude, py / magnitude, pz / magnitude)
        pressure = magnitude * u"Pa"
    end
    
    return AreaLoad(shells, pressure; 
        distribute_to=distribute_to, 
        interior_beams=interior_beams, 
        axis=axis, 
        direction=direction,
        id=id
    )
end

# Vector form for single shell
function AreaLoad(
    shell::ShellElement,
    load_vector::AbstractVector{<:Quantity};
    kwargs...
)
    AreaLoad([shell], load_vector; kwargs...)
end

"""
    nodal_forces(load::AreaLoad) -> Vector{Tuple{Node, Vector{Float64}}}

Compute equivalent nodal forces for an area load (FEM distribution).
Only used when distribute_to=:nodes.
"""
function nodal_forces(load::AreaLoad)
    load.distribute_to != :nodes && error("nodal_forces only valid for distribute_to=:nodes")
    
    p = ustrip(u"Pa", load.pressure)
    results = Tuple{Node, Vector{Float64}}[]
    
    for shell in load.shells
        area = shell.area
        total_force = p * area
        n_nodes = length(shell.nodes)
        force_per_node = total_force / n_nodes
        fvec = [force_per_node * d for d in load.direction]
        
        for node in shell.nodes
            push!(results, (node, fvec))
        end
    end
    
    return results
end

# =============================================================================
# Shell Self-Weight
# =============================================================================

"""
    SelfWeight(shells; g=9.81u"m/s^2")
    SelfWeight(shells, g)

Create an AreaLoad representing the self-weight of shell elements.

Computes pressure from shell properties: `p = ρ × t × g`

# Arguments
- `shells`: Shell elements (must have ρ defined)
- `g`: Gravitational acceleration (default: 9.81 m/s²)

# Example
```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
shells = Shell(corners, section)

# Self-weight load
sw = SelfWeight(shells)

# Or with custom g
sw = SelfWeight(shells; g=10u"m/s^2")
```
"""
function SelfWeight(
    shells::Vector{<:ShellElement};
    g::Quantity = 9.81u"m/s^2",
    id::Symbol = :selfweight
)
    isempty(shells) && error("SelfWeight requires at least one shell element")
    
    # Get thickness and density from first shell (assume uniform)
    first_shell = shells[1]
    t = first_shell.thickness  # meters (Float64)
    ρ = first_shell.ρ          # kg/m³ (Float64)
    
    # Compute self-weight pressure: p = ρ * t * g
    g_val = ustrip(u"m/s^2", g)
    pressure = ρ * t * g_val  # (kg/m³) * m * (m/s²) = kg/(m·s²) = Pa
    
    AreaLoad(shells, pressure * u"Pa"; direction=(0.0, 0.0, -1.0), id=id)
end

# Convenience: single shell
function SelfWeight(shell::ShellElement; kwargs...)
    SelfWeight([shell]; kwargs...)
end

# populate_load! for AreaLoad is defined in Model/preprocessing.jl
# to avoid circular dependency (loads.jl is included before model.jl)