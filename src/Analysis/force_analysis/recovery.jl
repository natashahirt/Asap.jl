# =============================================================================
# ElementInternalForces  —  struct
# =============================================================================

"""
    ElementInternalForces

Internal force results for a frame element (beam/column).

# Fields
- `element::Element` - Reference to the frame element
- `resolution::Integer` - Number of sampling points
- `x::Vector{Float64}` - Position along element [m]
- `P::Vector{Float64}` - Axial force [N]
- `My::Vector{Float64}` - Bending moment about local Y [N·m]
- `Vy::Vector{Float64}` - Shear force in local Y [N]
- `Mz::Vector{Float64}` - Bending moment about local Z [N·m]
- `Vz::Vector{Float64}` - Shear force in local Z [N]
"""
struct ElementInternalForces
    element::Element
    resolution::Integer
    x::Vector{Float64}
    P::Vector{Float64}
    My::Vector{Float64}
    Vy::Vector{Float64}
    Mz::Vector{Float64}
    Vz::Vector{Float64}
end

# =============================================================================
# ElementInternalForces  —  constructors
# =============================================================================

"""
    ElementInternalForces(element, model; resolution=20)

Internal force sampling for a frame element from a solved model.
"""
function ElementInternalForces(element::AbstractElement, model::Union{FrameModel, Model}; resolution = 20)
    
    dofs = etype2DOF[typeof(element)]
    ep = ElementProps(element)
    L = ep.L

    rng = range(0, L, resolution)

    # Use model.u directly (already Float64) + mul! to avoid intermediates
    u_e = model.u[element.globalID]
    Ku = Vector{Float64}(undef, 12)
    mul!(Ku, element.K, u_e)
    Flocal = Vector{Float64}(undef, 12)
    mul!(Flocal, element.R, Ku)
    @inbounds for i in 1:12
        Flocal[i] *= dofs[i]
    end

    Pstart  = -Flocal[1]
    Vystart =  Flocal[2]
    Mystart =  Flocal[6]
    Vzstart =  Flocal[3]
    Mzstart = -Flocal[5]

    P  = fill(Pstart, resolution)
    My = Vector{Float64}(undef, resolution)
    Vy = fill(Vystart, resolution)
    Mz = Vector{Float64}(undef, resolution)
    Vz = fill(Vzstart, resolution)
    @inbounds for i in 1:resolution
        My[i] = Vystart * rng[i] - Mystart
        Mz[i] = Vzstart * rng[i] - Mzstart
    end

    element_loads = get_elemental_loads(model)

    for load in element_loads[element.elementID]
        accumulate_force!(load, rng, P, My, Vy, Mz, Vz, ep)
    end

    return ElementInternalForces(element, resolution, collect(rng), P, My, Vy, Mz, Vz)
end

"""
    ElementInternalForces(element::Element, loads::Vector{<:ElementLoad}; resolution = 20)

Get internal force results for a given element from a set of loads.
Uses `mul!` with preallocated buffers to avoid R*K*u allocations.
"""
function ElementInternalForces(element::AbstractElement, loads::AbstractVector{<:AbstractLoad}; resolution = 20)
    
    dofs = etype2DOF[typeof(element)]
    ep = ElementProps(element)
    L = ep.L

    rng = range(0, L, resolution)

    # Fill uglobal directly (avoids 2 to_displacement_vec allocs + concat)
    _ug = Vector{Float64}(undef, 12)
    _fill_uglobal!(_ug, element)

    # Use mul! to avoid intermediate allocations for R * K * u
    Ku = Vector{Float64}(undef, 12)
    mul!(Ku, element.K, _ug)
    Flocal = Vector{Float64}(undef, 12)
    mul!(Flocal, element.R, Ku)
    @inbounds for i in 1:12
        Flocal[i] *= dofs[i]
    end

    Pstart  = -Flocal[1]
    Vystart =  Flocal[2]
    Mystart =  Flocal[6]
    Vzstart =  Flocal[3]
    Mzstart = -Flocal[5]

    P  = fill(Pstart, resolution)
    My = Vector{Float64}(undef, resolution)
    Vy = fill(Vystart, resolution)
    Mz = Vector{Float64}(undef, resolution)
    Vz = fill(Vzstart, resolution)
    @inbounds for i in 1:resolution
        My[i] = Vystart * rng[i] - Mystart
        Mz[i] = Vzstart * rng[i] - Mzstart
    end

    for load in loads
        accumulate_force!(load, rng, P, My, Vy, Mz, Vz, ep)
    end

    return ElementInternalForces(element, resolution, collect(rng), P, My, Vy, Mz, Vz)
end

