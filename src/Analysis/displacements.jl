"""
displacement function for the transverse translation and in-plane rotation for a GIVEN PLANE IN THE LCS OF AN ELEMENT:

u_xy= [ v₁
        θ₁
        v₂
        θ₂ ]

IE for the local XY plane:

v₁, v₂ are the start and end displacements in the local Y direction
θ₁, θ₂ are the start and end rotations in the **local Z** direction (ie rotation in plane of local XY)

Gives:

v_y(x) = N × u_xy (translational displacement in local Y at point x)

"""
function N(x::Float64, L::Float64)
    n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    n2 = x * (1 - x/L)^2
    n3 = 3(x/L)^2 - 2(x/L)^3
    n4 = x^2/L * (-1 + x/L)

    return [n1 n2 n3 n4]
end

"""
Axial displacement function: linear interpolation between start and end displacements
"""
function Naxial(x::Float64, L::Float64)
    n1 = 1 - x/L
    n2 = x / L

    return [n1 n2]
end

"""
    displacements(element::Element; n::Integer = 20)

Get the [3 × n] matrix where each column represents the local [x,y,z] displacement of the element from end forces.
Optimized to avoid broadcast+vcat allocation pattern.
"""
function unodal(element::Element; n::Integer = 20)

    # Fill uglobal directly (avoids 2 to_displacement_vec + concat allocs)
    _ug = Vector{Float64}(undef, 12)
    _fill_uglobal!(_ug, element)

    # In-place R * uglobal (avoids matvec alloc)
    _ul = Vector{Float64}(undef, 12)
    mul!(_ul, element.R, _ug)

    # Apply DOF mask in-place (avoids .* broadcast alloc)
    dofs = etype2DOF[typeof(element)]
    @inbounds for i in 1:12
        _ul[i] *= dofs[i]
    end

    L = ustrip(u"m", element.length)
    release_end_rotations!(_ul, element, L)

    uX1 = _ul[1];  uX2 = _ul[7]
    uY1 = _ul[2];  uY2 = _ul[6];  uY3 = _ul[8];  uY4 = _ul[12]
    uZ1 = _ul[3];  uZ2 = -_ul[5]; uZ3 = _ul[9];  uZ4 = -_ul[11]

    xrange = range(0, L, n)

    # Build displacement vectors directly (avoids n+1 small matrix allocations)
    dx = Vector{Float64}(undef, n)
    dy = Vector{Float64}(undef, n)
    dz = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        x = xrange[i]
        xL = x / L

        # Axial: linear interpolation
        dx[i] = (1 - xL) * uX1 + xL * uX2

        # Transverse (Hermite): N(x,L) * u
        n1 = 1 - 3xL^2 + 2xL^3
        n2 = x * (1 - xL)^2
        n3 = 3xL^2 - 2xL^3
        n4 = x^2/L * (-1 + xL)

        dy[i] = n1*uY1 + n2*uY2 + n3*uY3 + n4*uY4
        dz[i] = n1*uZ1 + n2*uZ2 + n3*uZ3 + n4*uZ4
    end

    dx, dy, dz
end

# =============================================================================
# accumulatedisp!  —  LineLoad
# =============================================================================

"""
Accumlate the internal forces cause by a line load to an element
"""
function accumulatedisp!(
    load::LineLoad, 
    xvals::AbstractVector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64},
    ep::ElementProps)

    v1 = ustrip(u"N/m", load.value[1])
    v2 = ustrip(u"N/m", load.value[2])
    v3 = ustrip(u"N/m", load.value[3])
    _, wy, wz = _rotate_to_local(ep.Rv, v1, v2, v3)

    @inbounds for i in eachindex(xvals)
        Dy[i] -= DLine(load.element, wy, ep.L, xvals[i], ep.E, ep.Ix)
        Dz[i] -= DLine(load.element, wz, ep.L, xvals[i], ep.E, ep.Iy)
    end
end


# =============================================================================
# accumulatedisp!  —  PointLoad
# =============================================================================

"""
Accumlate the internal forces cause by a point load to an element
"""
function accumulatedisp!(
    load::PointLoad, 
    xvals::AbstractVector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64},
    ep::ElementProps)

    frac = load.position

    v1 = ustrip(u"N", load.value[1])
    v2 = ustrip(u"N", load.value[2])
    v3 = ustrip(u"N", load.value[3])
    _, py, pz = _rotate_to_local(ep.Rv, v1, v2, v3)

    @inbounds for i in eachindex(xvals)
        Dy[i] -= DPoint(load.element, py, ep.L, xvals[i], frac, ep.E, ep.Ix)
        Dz[i] -= DPoint(load.element, pz, ep.L, xvals[i], frac, ep.E, ep.Iy)
    end
end


# =============================================================================
# accumulatedisp!  —  TributaryLoad
# =============================================================================

"""
Accumulate the displacement caused by a TributaryLoad to an element.
TributaryLoad represents a piecewise linear distributed load based on tributary widths.
The release-specific point-load deflection [`DPoint`](@ref) is integrated
against the load intensity, so the result is exact for every segment.
"""
function accumulatedisp!(
    load::TributaryLoad, 
    xvals::AbstractVector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64},
    ep::ElementProps)

    element = load.element
    L = ep.L

    d1, d2, d3 = load.direction
    _, uy, uz = _rotate_to_local(ep.Rv, d1, d2, d3)
    w = intensities(load)
    positions = load.positions

    for i in 1:(length(positions) - 1)
        a = positions[i] * L
        b = positions[i + 1] * L
        (b - a) < 1e-12 && continue
        wa, wb = w[i], w[i + 1]
        (abs(wa) + abs(wb)) < 1e-12 && continue
        wf(ξ) = wa + (wb - wa) * (ξ - a) / (b - a)
        @inbounds for j in eachindex(xvals)
            x = xvals[j]
            # DPoint is a polynomial in the load position on each side of x.
            fy(ξ) = wf(ξ) * DPoint(element, 1.0, L, x, ξ / L, ep.E, ep.Ix)
            fz(ξ) = wf(ξ) * DPoint(element, 1.0, L, x, ξ / L, ep.E, ep.Iy)
            lo, hi = a, b
            mid = clamp(x, lo, hi)
            Iy = _gauss3_integral(fy, lo, mid) + _gauss3_integral(fy, mid, hi)
            Iz = _gauss3_integral(fz, lo, mid) + _gauss3_integral(fz, mid, hi)
            Dy[j] -= uy * Iy
            Dz[j] -= uz * Iz
        end
    end
end

"""
    _gauss3_integral(f, lo, hi) -> Float64

Three-point Gauss–Legendre integral of `f` over `[lo, hi]`; exact for a
polynomial of degree ≤ 5.
"""
@inline function _gauss3_integral(f, lo, hi)
    hi - lo < 1e-15 && return 0.0
    h = (hi - lo) / 2
    m = (hi + lo) / 2
    s = 0.0
    @inbounds for k in 1:3
        s += _GAUSS3_W[k] * f(m + h * _GAUSS3_ξ[k])
    end
    return h * s
end


# =============================================================================
# accumulatedisp!  —  GravityLoad
# =============================================================================

"""
Accumulate the displacement cause by a GravityLoad (self-weight) to an element.
GravityLoad applies a distributed load equal to ρ * A * g in the global -Z direction.
"""
function accumulatedisp!(
    load::GravityLoad, 
    xvals::AbstractVector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64},
    ep::ElementProps)

    element = load.element

    ρ = ustrip(u"kg/m^3", element.section.ρ)
    A = ustrip(u"m^2", element.section.A)
    g = ustrip(u"m/s^2", load.factor)
    w_mag = ρ * A * g

    _, wy, wz = _rotate_to_local(ep.Rv, 0.0, 0.0, -w_mag)

    @inbounds for i in eachindex(xvals)
        Dy[i] -= DLine(element, wy, ep.L, xvals[i], ep.E, ep.Ix)
        Dz[i] -= DLine(element, wz, ep.L, xvals[i], ep.E, ep.Iy)
    end
end


# =============================================================================
# ElementDisplacements
# =============================================================================

"""
    ulocal(element::Element, model; resolution = 20)

Get the [3 × resolution] matrix of xyz displacements in LCS
"""
function ulocal(element::Element, model::Union{FrameModel, Model}; resolution = 20)
    ep = ElementProps(element)
    L = ep.L

    rng = range(0, L, resolution)

    Dx, Dy, Dz = unodal(element; n = resolution)

    element_loads = get_elemental_loads(model)

    for load in element_loads[element.elementID]
        accumulatedisp!(load, rng, Dy, Dz, ep)
    end

    return [Dx'; Dy'; Dz']
end

"""
    uglobal(element::Element, model; resolution = 20)

Get the [3 × resolution] matrix of xyz displacements in GCS
"""
function uglobal(element::Element, model::Union{FrameModel, Model}; resolution = 20)
    ep = ElementProps(element)
    L = ep.L

    rng = range(0, L, resolution)

    Dx, Dy, Dz = unodal(element; n = resolution)

    element_loads = get_elemental_loads(model)

    for load in element_loads[element.elementID]
        accumulatedisp!(load, rng, Dy, Dz, ep)
    end

    # Rotate local displacements to global using LCS vectors (avoids hcat + comprehension)
    lx, ly, lz = element.LCS[1], element.LCS[2], element.LCS[3]
    Dglobal = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        dx, dy, dz = Dx[j], Dy[j], Dz[j]
        Dglobal[1,j] = dx*lx[1] + dy*ly[1] + dz*lz[1]
        Dglobal[2,j] = dx*lx[2] + dy*ly[2] + dz*lz[2]
        Dglobal[3,j] = dx*lx[3] + dy*ly[3] + dz*lz[3]
    end
    return Dglobal
end

struct ElementDisplacements
    element::Element
    resolution::Integer
    x::Vector{Float64}
    ulocal::Matrix{Float64}
    uglobal::Matrix{Float64}
    basepositions::Matrix{Float64}
end

"""
    ElementDisplacements(element::Element, model; resolution = 20)

Get the local/global displacements of an element
"""
function ElementDisplacements(element::AbstractElement, model::Union{FrameModel, Model}; resolution = 20)
    element_loads = get_elemental_loads(model)
    return ElementDisplacements(element, element_loads[element.elementID]; resolution=resolution)
end

"""
    ElementDisplacements(element, loads; resolution=20)

Get displacements from pre-computed element loads (avoids rebuilding the loads map).
"""
function ElementDisplacements(element::AbstractElement, loads_for_elem::AbstractVector{<:AbstractLoad}; resolution = 20)
    ep = ElementProps(element)
    L = ep.L

    rng = range(0, L, resolution)

    Dx, Dy, Dz = unodal(element; n = resolution)

    for load in loads_for_elem
        accumulatedisp!(load, rng, Dy, Dz, ep)
    end

    # Build local displacement matrix directly from vectors (avoids [Dx'; Dy'; Dz'] allocation)
    D = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        D[1,j] = Dx[j]; D[2,j] = Dy[j]; D[3,j] = Dz[j]
    end

    # Rotate local displacements to global using LCS vectors
    lx, ly, lz = element.LCS[1], element.LCS[2], element.LCS[3]
    Dglobal = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        dx, dy, dz = Dx[j], Dy[j], Dz[j]
        Dglobal[1,j] = dx*lx[1] + dy*ly[1] + dz*lz[1]
        Dglobal[2,j] = dx*lx[2] + dy*ly[2] + dz*lz[2]
        Dglobal[3,j] = dx*lx[3] + dy*ly[3] + dz*lz[3]
    end

    # Base positions along element axis
    px = ustrip(u"m", element.nodeStart.position[1])
    py = ustrip(u"m", element.nodeStart.position[2])
    pz = ustrip(u"m", element.nodeStart.position[3])
    ax_x, ax_y, ax_z = element.LCS[1][1], element.LCS[1][2], element.LCS[1][3]
    xvec = collect(rng)
    basepoints = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        basepoints[1,j] = px + ax_x * rng[j]
        basepoints[2,j] = py + ax_y * rng[j]
        basepoints[3,j] = pz + ax_z * rng[j]
    end

    return ElementDisplacements(element, resolution, xvec, D, Dglobal, basepoints)
end

function ElementDisplacements(elements::AbstractVector{<:AbstractElement}, model::Union{FrameModel, Model}; resolution = 20)

    n_elem = length(elements)
    per_elem_res = max(Int(round(resolution / n_elem)), 2)
    total_pts = per_elem_res * n_elem

    xstore     = Vector{Float64}(undef, total_pts)
    ulocal_all = Matrix{Float64}(undef, 3, total_pts)
    uglobal_all = Matrix{Float64}(undef, 3, total_pts)
    basepoint_all = Matrix{Float64}(undef, 3, total_pts)

    element_loads = get_elemental_loads(model)
    k = 0

    for element in elements
        ep = ElementProps(element)
        L = ep.L

        rng = range(0, L, per_elem_res)

        Dx, Dy, Dz = unodal(element; n = per_elem_res)

        for load in element_loads[element.elementID]
            accumulatedisp!(load, rng, Dy, Dz, ep)
        end

        # Rotate local → global using LCS vectors
        lx, ly, lz = element.LCS[1], element.LCS[2], element.LCS[3]
        
        # Base position for this element
        px = ustrip(u"m", element.nodeStart.position[1])
        py = ustrip(u"m", element.nodeStart.position[2])
        pz = ustrip(u"m", element.nodeStart.position[3])
        ax_x, ax_y, ax_z = lx[1], lx[2], lx[3]  # local X axis = element axis

        x_offset = k > 0 ? xstore[k] : 0.0
        @inbounds for j in 1:per_elem_res
            k += 1
            xstore[k] = x_offset + rng[j]
            
            # Local displacements
            dx, dy, dz = Dx[j], Dy[j], Dz[j]
            ulocal_all[1,k] = dx
            ulocal_all[2,k] = dy
            ulocal_all[3,k] = dz
            
            # Global displacements (rotate LCS → GCS)
            uglobal_all[1,k] = dx*lx[1] + dy*ly[1] + dz*lz[1]
            uglobal_all[2,k] = dx*lx[2] + dy*ly[2] + dz*lz[2]
            uglobal_all[3,k] = dx*lx[3] + dy*ly[3] + dz*lz[3]
            
            # Base positions along element axis
            basepoint_all[1,k] = px + ax_x * rng[j]
            basepoint_all[2,k] = py + ax_y * rng[j]
            basepoint_all[3,k] = pz + ax_z * rng[j]
        end
    end

    resize!(xstore, k)
    ulocal_out = ulocal_all[:, 1:k]
    uglobal_out = uglobal_all[:, 1:k]
    basepoint_out = basepoint_all[:, 1:k]
    return ElementDisplacements(elements[1], resolution, xstore, ulocal_out, uglobal_out, basepoint_out)
end

"""
    displacements(model, increment)

Get the displacements of all elements in a model.

# Arguments
- `model` - Structural model (FrameModel or Model)
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)
"""
function displacements(model::FrameModel, increment)
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)

    results = Vector{ElementDisplacements}()
    sizehint!(results, length(model.elements))
    for element in model.elements
        L = ustrip(u"m", element.length)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, ElementDisplacements(element, model; resolution = n))
    end

    return results
end

function displacements(model::Model, increment)
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)
    
    if isempty(model.frame_elements)
        return ElementDisplacements[]
    end
    
    results = Vector{ElementDisplacements}()
    ids = groupbyid(model.frame_elements)
    sizehint!(results, length(ids))

    for id in ids
        elements = model.frame_elements[id]
        L = sum(ustrip(u"m", el.length) for el in elements)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, ElementDisplacements(elements, model; resolution = n))
    end

    return results
end



# =============================================================================
# Vertical deflection from an arbitrary displacement vector and load case
# =============================================================================

"""
    deflection_axis_local(element) -> Union{Nothing, NTuple{2, Float64}}

Components `(n_y, n_z)`, along the element's local y and z axes, of the unit
vector `n` that is perpendicular to the element axis and lies in the vertical
plane containing it:

    n = (ẑ − (ẑ·x̂) x̂) / |ẑ − (ẑ·x̂) x̂|

with `ẑ` the global +Z axis and `x̂` the element axis. `n` depends only on the
element's direction, not on its roll angle Ψ, so a horizontal beam returns the
global vertical for the default Ψ = π/2 (local y vertical), for Ψ = 0 (local z
vertical), and for any other roll. Returns `nothing` for a vertical element,
which has no such plane.
"""
function deflection_axis_local(element::AbstractElement)
    lx, ly, lz = element.LCS[1], element.LCS[2], element.LCS[3]
    c = lx[3]
    nx = -c * lx[1]
    ny = -c * lx[2]
    nz = 1.0 - c * lx[3]
    nn = sqrt(nx^2 + ny^2 + nz^2)
    nn < 1e-6 && return nothing
    n_y = (nx * ly[1] + ny * ly[2] + nz * ly[3]) / nn
    n_z = (nx * lz[1] + ny * lz[2] + nz * lz[3]) / nn
    return (n_y, n_z)
end

"""
    element_deflection_profile(element, u, loads; resolution=20) -> Vector{Float64}

Transverse deflection (m) of `element` at `resolution` evenly spaced stations
from the start node to the end node, for the global displacement vector `u` and
the member loads `loads` that produced it.

Each station is the sum of
1. the Hermite cubic through the element's end translations and slopes, read
   from `u` after the `etype2DOF` mask and [`release_end_rotations!`](@ref), and
2. the member-load particular solution of every load in `loads`
   ([`accumulatedisp!`](@ref)), which is zero at both ends.

Without part 2 a span modelled as one element reports no deflection from its own
load. `loads` must be the element loads of the same case as `u`, all on
`element`; a load on another element throws.

The returned value is the component along [`deflection_axis_local`](@ref): the
deflection in the vertical plane through the member, positive toward global +Z,
which is the deflection that span/δ serviceability limits measure. It is the
absolute displacement, so it includes movement of the member's end nodes. A
vertical element has no vertical plane; it returns the transverse magnitude
`√(Dy² + Dz²)` instead.
"""
function element_deflection_profile(element::AbstractElement,
                                    u::AbstractVector{<:Real},
                                    loads::AbstractVector{<:AbstractLoad};
                                    resolution::Int = 20)
    resolution >= 2 || throw(ArgumentError(
        "element_deflection_profile: resolution must be at least 2 (got $resolution)"))
    ep = ElementProps(element)
    L = ep.L
    L > 0 || throw(ArgumentError(
        "element_deflection_profile: element $(element.elementID) has zero length"))

    gid = element.globalID
    ug = Vector{Float64}(undef, 12)
    @inbounds for i in 1:12
        ug[i] = u[gid[i]]
    end
    ul = Vector{Float64}(undef, 12)
    mul!(ul, element.R, ug)
    dofs = etype2DOF[typeof(element)]
    @inbounds for i in 1:12
        ul[i] *= dofs[i]
    end
    release_end_rotations!(ul, element, L)

    # Local xy plane: (v, θz); local xz plane: (w, −θy). See `unodal`.
    uY1 = ul[2]; uY2 = ul[6];  uY3 = ul[8]; uY4 = ul[12]
    uZ1 = ul[3]; uZ2 = -ul[5]; uZ3 = ul[9]; uZ4 = -ul[11]

    xs = range(0.0, L, resolution)
    Dy = Vector{Float64}(undef, resolution)
    Dz = Vector{Float64}(undef, resolution)
    @inbounds for i in 1:resolution
        x = xs[i]
        xL = x / L
        n1 = 1 - 3xL^2 + 2xL^3
        n2 = x * (1 - xL)^2
        n3 = 3xL^2 - 2xL^3
        n4 = x^2 / L * (-1 + xL)
        Dy[i] = n1 * uY1 + n2 * uY2 + n3 * uY3 + n4 * uY4
        Dz[i] = n1 * uZ1 + n2 * uZ2 + n3 * uZ3 + n4 * uZ4
    end

    for load in loads
        load.element === element || throw(ArgumentError(
            "element_deflection_profile: a load on element $(load.element.elementID) " *
            "was passed for element $(element.elementID)"))
        accumulatedisp!(load, xs, Dy, Dz, ep)
    end

    axis = deflection_axis_local(element)
    out = Vector{Float64}(undef, resolution)
    if axis === nothing
        @inbounds for i in 1:resolution
            out[i] = hypot(Dy[i], Dz[i])
        end
    else
        n_y, n_z = axis
        @inbounds for i in 1:resolution
            out[i] = n_y * Dy[i] + n_z * Dz[i]
        end
    end
    return out
end

"""
    _area_load_tributaries!(load::AreaLoad) -> Vector{TributaryLoad}

Member loads an `AreaLoad` with `distribute_to` frame elements puts on those
elements. Same conversion and cache as `populate_load!(model, ::AreaLoad)`; the
cached loads follow `load.pressure` through `AreaLoad`'s `setproperty!`.
"""
function _area_load_tributaries!(load::AreaLoad)
    if isnothing(load._tributary_loads)
        axis_vec = isnothing(load.axis) ? nothing : [load.axis[1], load.axis[2]]
        load._tributary_loads = _shell_to_tributary_loads(
            load.shells, load.distribute_to, load.pressure;
            axis = axis_vec, direction = load.direction,
            interior_beams = load.interior_beams)
    end
    return load._tributary_loads
end

"""
    frame_loads_by_element(model, loads) -> Vector{Vector{AbstractLoad}}

Member loads in `loads`, grouped by frame element in model order
(`model.elements` for a `FrameModel`, `model.frame_elements` for a `Model`).
Loads are matched to elements by identity, not `elementID`. An `AreaLoad`
distributed to frame elements contributes its tributary member loads; nodal
loads, nodal area loads, and shell-element loads carry no member load and are
skipped. A load on a frame element that is not in the model throws.
"""
function frame_loads_by_element(model::Union{FrameModel, Model},
                                loads::AbstractVector{<:AbstractLoad})
    elements = model isa FrameModel ? model.elements : model.frame_elements
    index = IdDict{Any, Int}()
    for (i, el) in enumerate(elements)
        index[el] = i
    end
    out = [AbstractLoad[] for _ in elements]
    function add!(load)
        i = get(index, load.element, 0)
        i == 0 && throw(ArgumentError(
            "frame_loads_by_element: a $(nameof(typeof(load))) is on a frame element " *
            "that is not in the model"))
        push!(out[i], load)
    end
    for load in loads
        if load isa AreaLoad
            load.distribute_to isa Symbol && continue
            foreach(add!, _area_load_tributaries!(load))
        elseif hasproperty(load, :element) && load.element isa FrameElement
            add!(load)
        end
    end
    return out
end

"""
    element_deflection_profiles(model, u, loads; resolution=20) -> Vector{Vector{Float64}}

[`element_deflection_profile`](@ref) for every frame element, in model order,
for the displacement vector `u` of one load case and that case's `loads`.
Profiles of separate cases on the same stiffness add station by station.
"""
function element_deflection_profiles(model::Union{FrameModel, Model},
                                     u::AbstractVector{<:Real},
                                     loads::AbstractVector{<:AbstractLoad};
                                     resolution::Int = 20)
    length(u) == model.nDOFs || throw(DimensionMismatch(
        "element_deflection_profiles: u has length $(length(u)), the model has " *
        "$(model.nDOFs) DOFs"))
    elements = model isa FrameModel ? model.elements : model.frame_elements
    by_el = frame_loads_by_element(model, loads)
    return [element_deflection_profile(el, u, by_el[i]; resolution = resolution)
            for (i, el) in enumerate(elements)]
end

"""
    element_max_deflections(model, u, loads; resolution=20) -> Dict{Int, Float64}

Largest absolute [`element_deflection_profile`](@ref) value (m) of each frame
element, keyed by `elementID` (the largest over elements that share an ID), for
the displacement vector `u` of one load case and that case's `loads`.

`u` can be any solved vector (service dead, live, or their sum) and the model is
not mutated, apart from caching the tributary loads of an `AreaLoad`. Pass the
loads of the same case: the member-load term is what gives a span modelled as
one element its own deflection. For the deflection of a sum of cases, add the
profiles from [`element_deflection_profiles`](@ref) before taking the maximum.
"""
function element_max_deflections(model::Union{FrameModel, Model},
                                 u::AbstractVector{<:Real},
                                 loads::AbstractVector{<:AbstractLoad};
                                 resolution::Int = 20)
    elements = model isa FrameModel ? model.elements : model.frame_elements
    profiles = element_deflection_profiles(model, u, loads; resolution = resolution)
    result = Dict{Int, Float64}()
    sizehint!(result, length(elements))
    for (i, el) in enumerate(elements)
        δ = maximum(abs, profiles[i])
        eid = el.elementID
        result[eid] = max(get(result, eid, 0.0), δ)
    end
    return result
end
