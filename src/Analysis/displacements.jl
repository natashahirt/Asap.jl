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
Each segment is approximated as a point load at its centroid.
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
    pressure = load._pressure_Pa
    widths_m = load._widths_m
    positions = load.positions
    n_seg = length(positions) - 1

    for i in 1:n_seg
        s1, s2 = positions[i], positions[i+1]
        w1, w2 = widths_m[i], widths_m[i+1]

        (s2 - s1) < 1e-12 && continue
        (w1 + w2) < 1e-12 && continue

        a = s1 * L
        b = s2 * L
        seg_len = b - a
        
        w_avg = pressure * (w1 + w2) / 2
        F_total = w_avg * seg_len
        frac = (a + b) / 2 / L

        g1, g2, g3 = d1 * F_total, d2 * F_total, d3 * F_total
        _, py, pz = _rotate_to_local(ep.Rv, g1, g2, g3)

        @inbounds for j in eachindex(xvals)
            Dy[j] -= DPoint(element, py, L, xvals[j], frac, ep.E, ep.Ix)
            Dz[j] -= DPoint(element, pz, L, xvals[j], frac, ep.E, ep.Iy)
        end
    end
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
# element_max_deflections  —  nodal-only max |δ| per element from arbitrary u
# =============================================================================

"""
    element_max_deflections(model, u; resolution=20)

Return a `Dict{Int, Float64}` mapping each frame element's `elementID` to
its maximum absolute transverse deflection (meters) along the element,
computed from the displacement vector `u` via Hermite cubic interpolation.

Unlike `ElementDisplacements`, this function reads directly from `u` rather
than from `node.displacement`, so it can be used with any displacement vector
(e.g. service-load D-only, L-only, or D+L) without mutating the model.

Only nodal shape-function interpolation is used (no load-based corrections),
which is accurate when spans are meshed with multiple elements.

The transverse deflection is measured in the element's **local Z** direction
(strong-axis bending plane for standard floor beams).

# Arguments
- `model` — `FrameModel` or `Model`
- `u::Vector{Float64}` — global displacement vector (length = nDOFs)
- `resolution::Int` — sampling points per element (default 20)
"""
function element_max_deflections(
    model::Union{FrameModel, Model},
    u::Vector{Float64};
    resolution::Int = 20,
)
    elements = model isa FrameModel ? model.elements : model.frame_elements
    result = Dict{Int, Float64}()
    sizehint!(result, length(elements))

    _ug = Vector{Float64}(undef, 12)
    _ul = Vector{Float64}(undef, 12)

    for el in elements
        gid = el.globalID
        @inbounds for i in 1:12
            _ug[i] = u[gid[i]]
        end
        mul!(_ul, el.R, _ug)

        dofs = etype2DOF[typeof(el)]
        @inbounds for i in 1:12
            _ul[i] *= dofs[i]
        end

        L = ustrip(u"m", el.length)
        L <= 0 && continue

        uZ1 = _ul[3]; uZ2 = -_ul[5]; uZ3 = _ul[9]; uZ4 = -_ul[11]

        δ_max = 0.0
        @inbounds for i in 0:(resolution - 1)
            xL = i / (resolution - 1)
            x  = xL * L
            n1 = 1 - 3xL^2 + 2xL^3
            n2 = x * (1 - xL)^2
            n3 = 3xL^2 - 2xL^3
            n4 = x^2/L * (-1 + xL)
            dz = abs(n1*uZ1 + n2*uZ2 + n3*uZ3 + n4*uZ4)
            δ_max = max(δ_max, dz)
        end

        eid = el.elementID
        prev = get(result, eid, 0.0)
        result[eid] = max(prev, δ_max)
    end

    return result
end
