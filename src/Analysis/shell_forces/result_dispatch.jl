# ============================================================================
# Unified InternalForces Interface
# ============================================================================

"""
    InternalForces(element, model; kwargs...)

Unified interface for internal force recovery. Dispatches to:
- `ElementInternalForces` for frame elements (beams, columns, trusses)
- `ShellInternalForces` for shell elements

# Example
```julia
model = Model(nodes, frames, shells, [load])
process!(model)
solve!(model)

beam_forces = InternalForces(frames[1], model)   # → ElementInternalForces
shell_forces = InternalForces(shells[1], model)  # → ShellInternalForces
```
"""
function InternalForces(elem::FrameElement, model::AbstractModel; resolution=20)
    ElementInternalForces(elem, model; resolution=resolution)
end

function InternalForces(elem::ShellElement, model::AbstractModel)
    ShellInternalForces(elem, model.u)
end

"""
    InternalForces(model; increment=0.1)

Get internal forces for all elements in a model.

# Returns
Named tuple: `(frames=Vector{ElementInternalForces}, shells=Vector{ShellInternalForces})`
"""
function InternalForces(model::Model; increment=0.1)
    frame_forces = ElementInternalForces(model, increment)
    shell_forces = ShellInternalForces(model)
    (frames=frame_forces, shells=shell_forces)
end

function InternalForces(model::FrameModel; increment=0.1)
    (frames=ElementInternalForces(model, increment), shells=ShellInternalForces[])
end

function InternalForces(model::ShellModel)
    (frames=ElementInternalForces[], shells=ShellInternalForces(model))
end

# ============================================================================
# Unified Displacements Interface
# ============================================================================

"""
    Displacements(element, model)

Get displacement for a single element.

# Returns
- `ShellDisplacements` for shell elements
- For frame elements, use `displacements(element, model)` instead

# Example
```julia
disp = Displacements(shell1, model)
println("w = \$(disp.w) m")
```
"""
function Displacements(elem::ShellElement, model::AbstractModel)
    ShellDisplacements(elem, model.u)
end

"""
    Displacements(model)

Get displacements for all elements in a model.

# Returns
Named tuple: `(shells=Vector{ShellDisplacements},)`

# Example
```julia
disps = Displacements(model)
disps.shells  # All shell centroid displacements
```
"""
function Displacements(model::Model)
    (shells=ShellDisplacements(model),)
end

function Displacements(model::ShellModel)
    (shells=ShellDisplacements(model),)
end
