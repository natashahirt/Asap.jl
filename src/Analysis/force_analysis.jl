
# =============================================================================
# Helper: rotate a 3-vector from global to local with beam sign convention
# =============================================================================

"""
    _rotate_to_local(Rv, g1, g2, g3)

Rotate global vector (g1,g2,g3) into the local frame using the 3×3 rotation
sub-block `Rv = @view element.R[1:3,1:3]`.  Returns `(wx, wy, wz)` with the
beam sign convention (y,z flipped).
"""
@inline function _rotate_to_local(Rv, g1, g2, g3)
    @inbounds begin
        wx =  Rv[1,1]*g1 + Rv[1,2]*g2 + Rv[1,3]*g3
        wy = -(Rv[2,1]*g1 + Rv[2,2]*g2 + Rv[2,3]*g3)
        wz = -(Rv[3,1]*g1 + Rv[3,2]*g2 + Rv[3,3]*g3)
    end
    return wx, wy, wz
end

# =============================================================================
# ElementProps — pre-stripped element properties for hot-path loops
# =============================================================================

"""
    ElementProps

Pre-stripped (unitless Float64) element properties, computed once per element
and passed into `accumulate_force!` / `accumulatedisp!` to avoid repeated
`ustrip` calls inside inner loops.
"""
struct ElementProps
    Rv::SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}
    L::Float64
    E::Float64
    Ix::Float64   # strong-axis I (I₂ / Iy in Asap convention)
    Iy::Float64   # weak-axis I  (I₃ / Iz in Asap convention)
end

"""Pre-strip element properties once for hot-path force/displacement loops."""
@inline function ElementProps(element::AbstractElement)
    ElementProps(
        (@view element.R[1:3, 1:3]),
        ustrip(u"m", element.length),
        ustrip(u"Pa", element.section.E),
        ustrip(u"m^4", element.section.Ix),
        ustrip(u"m^4", element.section.Iy),
    )
end

# =============================================================================
# _fill_uglobal! — fill 12-vec with nodal displacements (avoids to_displacement_vec + concat)
# =============================================================================

"""
    _fill_uglobal!(buf, element)

Fill a 12-element buffer with nodal displacements directly from element node
data.  DOFs 1-3 are translations (m), DOFs 4-6 are rotations (rad).
Avoids two `to_displacement_vec` allocations + `[start; end]` concatenation.
"""
@inline function _fill_uglobal!(buf::Vector{Float64}, element::AbstractElement)
    ds = element.nodeStart.displacement
    de = element.nodeEnd.displacement
    @inbounds begin
        buf[1]  = ustrip(u"m", ds[1])
        buf[2]  = ustrip(u"m", ds[2])
        buf[3]  = ustrip(u"m", ds[3])
        buf[4]  = ustrip(ds[4])
        buf[5]  = ustrip(ds[5])
        buf[6]  = ustrip(ds[6])
        buf[7]  = ustrip(u"m", de[1])
        buf[8]  = ustrip(u"m", de[2])
        buf[9]  = ustrip(u"m", de[3])
        buf[10] = ustrip(de[4])
        buf[11] = ustrip(de[5])
        buf[12] = ustrip(de[6])
    end
    return buf
end

# =============================================================================
# accumulate_force!  —  TributaryLoad
# =============================================================================

"""
Accumulate the internal forces caused by a TributaryLoad to the current element.
TributaryLoad represents a piecewise linear distributed load based on tributary widths.
Each segment between breakpoints is approximated as a partial uniform load.
"""
function accumulate_force!(load::TributaryLoad, 
    xvals::AbstractVector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64},
    ep::ElementProps)

    element = load.element
    L = ep.L
    Rv = ep.Rv
    
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
        
        w_avg = pressure * (w1 + w2) / 2
        g1, g2, g3 = d1 * w_avg, d2 * w_avg, d3 * w_avg
        wx, wy, wz = _rotate_to_local(Rv, g1, g2, g3)
        
        @inbounds for (j, x) in enumerate(xvals)
            P[j]  += _PLine_partial(wx, L, x, a, b)
            My[j] += _MLine_partial(element, wy, L, x, a, b)
            Vy[j] += _VLine_partial(element, wy, L, x, a, b)
            Mz[j] += _MLine_partial(element, wz, L, x, a, b)
            Vz[j] += _VLine_partial(element, wz, L, x, a, b)
        end
    end
end


# Helper: axial force from partial uniform load on [a, b]
@inline function _PLine_partial(w, L, x, a, b)
    F = w * (b - a)
    R_start = F * (L - (a + b) / 2) / L
    
    if x < a
        -R_start
    elseif x < b
        -(R_start + w * (x - a))
    else
        -(R_start + F)
    end
end

# =============================================================================
# End Condition Behavior Traits
# =============================================================================

abstract type EndConditionBehavior end
struct SimplySupportedBehavior <: EndConditionBehavior end
struct FixedFixedBehavior <: EndConditionBehavior end

_end_behavior(::Element{FreeFree}) = SimplySupportedBehavior()
_end_behavior(::Element{Joist}) = SimplySupportedBehavior()
_end_behavior(::Element{FixedFixed}) = FixedFixedBehavior()
_end_behavior(::Element{FreeFixed}) = SimplySupportedBehavior()
_end_behavior(::Element{FixedFree}) = SimplySupportedBehavior()

function _MLine_partial(elem::Element, w, L, x, a, b)
    _MLine_partial(_end_behavior(elem), w, L, x, a, b)
end

_MLine_partial(::SimplySupportedBehavior, w, L, x, a, b) = _MLine_partial_ss(w, L, x, a, b)
_MLine_partial(::FixedFixedBehavior, w, L, x, a, b) = _MLine_partial_ff(w, L, x, a, b)

@inline function _MLine_partial_ss(w, L, x, a, b)
    c = (a + b) / 2
    F = w * (b - a)
    Ra = F * (L - c) / L
    
    if x <= a
        Ra * x
    elseif x <= b
        Ra * x - w * (x - a)^2 / 2
    else
        Ra * x - F * (x - c)
    end
end

@inline function _MLine_partial_ff(w, L, x, a, b)
    seg_len = b - a
    c = (a + b) / 2
    F = w * seg_len
    
    Ma = F * c * (L - c)^2 / L^2 - F * seg_len^2 / 12 * (L - 2*c) / L
    Mb = F * (L - c) * c^2 / L^2 + F * seg_len^2 / 12 * (L - 2*c) / L
    
    Ra = (F * (L - c) + Ma - Mb) / L
    
    if x <= a
        Ra * x - Ma
    elseif x <= b
        Ra * x - Ma - w * (x - a)^2 / 2
    else
        Ra * x - Ma - F * (x - c)
    end
end

_MLine_partial_pf(w, L, x, a, b) = _MLine_partial_ss(w, L, x, a, b)
_MLine_partial_fp(w, L, x, a, b) = _MLine_partial_ss(w, L, x, a, b)

@inline function _VLine_partial(element, w, L, x, a, b)
    c = (a + b) / 2
    F = w * (b - a)
    Ra = F * (L - c) / L
    
    if x <= a
        Ra
    elseif x <= b
        Ra - w * (x - a)
    else
        Ra - F
    end
end

# =============================================================================
# accumulate_force!  —  LineLoad
# =============================================================================

"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulate_force!(load::LineLoad, 
    xvals::AbstractVector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64},
    ep::ElementProps)

    L = ep.L

    v1 = ustrip(u"N/m", load.value[1])
    v2 = ustrip(u"N/m", load.value[2])
    v3 = ustrip(u"N/m", load.value[3])
    wx, wy, wz = _rotate_to_local(ep.Rv, v1, v2, v3)

    @inbounds for (j, x) in enumerate(xvals)
        P[j]  += PLine(wx, L, x)
        My[j] += MLine(load.element, wy, L, x)
        Vy[j] += VLine(load.element, wy, L, x)
        Mz[j] += MLine(load.element, wz, L, x)
        Vz[j] += VLine(load.element, wz, L, x)
    end
end


# =============================================================================
# accumulate_force!  —  PointLoad
# =============================================================================

"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulate_force!(load::PointLoad, 
    xvals::AbstractVector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64},
    ep::ElementProps)

    L = ep.L
    frac = load.position

    v1 = ustrip(u"N", load.value[1])
    v2 = ustrip(u"N", load.value[2])
    v3 = ustrip(u"N", load.value[3])
    px, py, pz = _rotate_to_local(ep.Rv, v1, v2, v3)

    @inbounds for (j, x) in enumerate(xvals)
        P[j]  += PPoint(px, L, x, frac)
        My[j] += MPoint(load.element, py, L, x, frac)
        Vy[j] += VPoint(load.element, py, L, x, frac)
        Mz[j] += MPoint(load.element, pz, L, x, frac)
        Vz[j] += VPoint(load.element, pz, L, x, frac)
    end
end


# =============================================================================
# accumulate_force!  —  GravityLoad
# =============================================================================

"""
Accumulate the internal forces caused by a GravityLoad (self-weight) to the current element.
GravityLoad applies a distributed load equal to ρ * A * g in the global -Z direction.
"""
function accumulate_force!(load::GravityLoad, 
    xvals::AbstractVector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64},
    ep::ElementProps)

    element = load.element
    L = ep.L

    ρ = ustrip(u"kg/m^3", element.section.ρ)
    A = ustrip(u"m^2", element.section.A)
    g = ustrip(u"m/s^2", load.factor)
    w_mag = ρ * A * g

    wx, wy, wz = _rotate_to_local(ep.Rv, 0.0, 0.0, -w_mag)

    @inbounds for (j, x) in enumerate(xvals)
        P[j]  += PLine(wx, L, x)
        My[j] += MLine(element, wy, L, x)
        Vy[j] += VLine(element, wy, L, x)
        Mz[j] += MLine(element, wz, L, x)
        Vz[j] += VLine(element, wz, L, x)
    end
end


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

# =============================================================================
# ElementInternalForces — sizehinted batch constructors
# =============================================================================


"""
    get_elemental_loads(model) -> Vector{Vector{AbstractLoad}}

Return a per-element load map.  The result is lazily cached on the model
(`model._elemental_loads`) and invalidated automatically by `process!`.
"""
function get_elemental_loads(model::FrameModel)
    cache = model._elemental_loads
    cache !== nothing && return cache
    element_to_loads = [AbstractLoad[] for _ in 1:model.nElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        push!(element_to_loads[load.element.elementID], load)
    end
    model._elemental_loads = element_to_loads
    return element_to_loads
end

function get_elemental_loads(model::Model)
    cache = model._elemental_loads
    cache !== nothing && return cache
    element_to_loads = [AbstractLoad[] for _ in 1:model.nFrameElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        if load.element isa FrameElement
            push!(element_to_loads[load.element.elementID], load)
        end
    end
    model._elemental_loads = element_to_loads
    return element_to_loads
end

"""
    get_elemental_loads(loads, n_elements) -> Vector{Vector{AbstractLoad}}

Build a per-element load map from an arbitrary load vector (no model required).
Useful for post-processing displacement vectors from `solve(model, cases)` where
each case has its own load set.
"""
function get_elemental_loads(loads::Vector{<:AbstractLoad}, n_elements::Int)
    element_to_loads = [AbstractLoad[] for _ in 1:n_elements]
    for load in loads
        hasproperty(load, :element) || continue
        eid = load.element.elementID
        (1 <= eid <= n_elements) && push!(element_to_loads[eid], load)
    end
    return element_to_loads
end

"""Invalidate the cached elemental loads map (call when loads or elements change)."""
function clear_elemental_loads!(model::AbstractModel)
    hasproperty(model, :_elemental_loads) && (model._elemental_loads = nothing)
end

"""
    ElementInternalForces(elements::Vector{<:FrameElement}, model; resolution=20)

Get internal force results for a group of ordered elements that form a single physical element.
"""
function ElementInternalForces(elements::AbstractVector{<:AbstractElement}, model::Union{FrameModel, Model}; resolution = 20)
    
    per_elem_res = max(Int(round(resolution / length(elements))), 2)
    total_pts = per_elem_res * length(elements)

    xstore  = Vector{Float64}(undef, total_pts)
    pstore  = Vector{Float64}(undef, total_pts)
    mystore = Vector{Float64}(undef, total_pts)
    vystore = Vector{Float64}(undef, total_pts)
    mzstore = Vector{Float64}(undef, total_pts)
    vzstore = Vector{Float64}(undef, total_pts)

    # Scratch vectors reused per element
    P_tmp  = Vector{Float64}(undef, per_elem_res)
    My_tmp = Vector{Float64}(undef, per_elem_res)
    Vy_tmp = Vector{Float64}(undef, per_elem_res)
    Mz_tmp = Vector{Float64}(undef, per_elem_res)
    Vz_tmp = Vector{Float64}(undef, per_elem_res)
    xinc   = Vector{Float64}(undef, per_elem_res)

    element_loads_map = get_elemental_loads(model)
    k = 0

    # Preallocated mul! buffers for R*K*u
    _Ku_buf = Vector{Float64}(undef, 12)
    _Fl_buf = Vector{Float64}(undef, 12)

    for element in elements
        loadids = element_loads_map[element.elementID]

        dofs = etype2DOF[typeof(element)]
        ep = ElementProps(element)
        L = ep.L

        # Fill xinc
        rng = range(0, L, per_elem_res)
        @inbounds for j in 1:per_elem_res
            xinc[j] = rng[j]
        end

        # Use model.u directly (Float64) + mul! to avoid intermediates
        u_e = model.u[element.globalID]
        mul!(_Ku_buf, element.K, u_e)
        mul!(_Fl_buf, element.R, _Ku_buf)
        @inbounds for i in 1:12
            _Fl_buf[i] *= dofs[i]
        end

        Pstart  = -_Fl_buf[1]
        Vystart =  _Fl_buf[2]
        Mystart =  _Fl_buf[6]
        Vzstart =  _Fl_buf[3]
        Mzstart = -_Fl_buf[5]

        fill!(P_tmp, Pstart)
        fill!(Vy_tmp, Vystart)
        fill!(Vz_tmp, Vzstart)
        @inbounds for j in 1:per_elem_res
            My_tmp[j] = Vystart * xinc[j] - Mystart
            Mz_tmp[j] = Vzstart * xinc[j] - Mzstart
        end

        for load in loadids
            accumulate_force!(load, xinc, P_tmp, My_tmp, Vy_tmp, Mz_tmp, Vz_tmp, ep)
        end

        x_offset = k > 0 ? xstore[k] : 0.0
        @inbounds for j in 1:per_elem_res
            k += 1
            xstore[k]  = x_offset + xinc[j]
            pstore[k]  = P_tmp[j]
            mystore[k] = My_tmp[j]
            vystore[k] = Vy_tmp[j]
            mzstore[k] = Mz_tmp[j]
            vzstore[k] = Vz_tmp[j]
        end
    end

    resize!(xstore, k); resize!(pstore, k); resize!(mystore, k)
    resize!(vystore, k); resize!(mzstore, k); resize!(vzstore, k)

    return ElementInternalForces(elements[1], resolution, xstore, pstore, mystore, vystore, mzstore, vzstore)
end



"""
    ElementInternalForces(model::FrameModel, increment)

Get the internal forces of all frame elements in a FrameModel.

# Arguments
- `model` - FrameModel (frame elements only)
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)
"""
function ElementInternalForces(model::FrameModel, increment)
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)
    
    element_loads = get_elemental_loads(model)

    IF = [ElementInternalForces(e, model.u[e.globalID], element_loads[e.elementID], inc_m, etype2DOF[typeof(e)]) for e in model.elements]

    # Preallocate concatenated arrays (avoids vcat + splat allocations)
    total_pts = sum(length(IF[i].x) for i in eachindex(IF))
    L_total = sum(IF[i].L for i in eachindex(IF))
    x  = Vector{Float64}(undef, total_pts)
    P  = Vector{Float64}(undef, total_pts)
    My = Vector{Float64}(undef, total_pts)
    Vy = Vector{Float64}(undef, total_pts)
    Mz = Vector{Float64}(undef, total_pts)
    Vz = Vector{Float64}(undef, total_pts)

    offset = 0
    x_cum = 0.0
    @inbounds for i in eachindex(IF)
        eif = IF[i]
        n = length(eif.x)
        x_off = x_cum + eif.L / 2
        for j in 1:n
            k = offset + j
            x[k]  = eif.x[j] + x_off
            P[k]  = eif.P[j]
            My[k] = eif.My[j]
            Vy[k] = eif.Vy[j]
            Mz[k] = eif.Mz[j]
            Vz[k] = eif.Vz[j]
        end
        offset += n
        x_cum += eif.L
    end

    return ElementInternalForces(nothing, L_total, x, P, My, Vy, Mz, Vz)
end

"""
    ElementInternalForces(model::Model, increment)

Get the internal forces of all frame elements in a Model.

# Arguments
- `model` - Model (may contain both frame and shell elements)
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)

# Returns
Vector of `ElementInternalForces`, one per frame element.
"""
function ElementInternalForces(model::Model, increment)
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)

    if isempty(model.frame_elements)
        return ElementInternalForces[]
    end

    results = Vector{ElementInternalForces}()
    sizehint!(results, length(model.frame_elements))
    element_loads_map = get_elemental_loads(model)

    for element in model.frame_elements
        L = ustrip(u"m", element.length)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, ElementInternalForces(element, element_loads_map[element.elementID]; resolution = n))
    end

    return results
end

struct ForceEnvelopes
    element::Element
    resolution::Integer
    x::Vector{Float64}
    Plow::Vector{Float64}
    Phigh::Vector{Float64}
    Mylow::Vector{Float64}
    Myhigh::Vector{Float64}
    Vylow::Vector{Float64}
    Vyhigh::Vector{Float64}
    Mzlow::Vector{Float64}
    Mzhigh::Vector{Float64}
    Vzlow::Vector{Float64}
    Vzhigh::Vector{Float64}
end

"""
    load_envelopes(model, loads::Vector{Vector{<:AbstractLoad}})

Get the high/low internal forces for a series of external loads
"""
function load_envelopes(model::Union{FrameModel, Model}, loads::Vector{Vector{<:AbstractLoad}}, increment::Real)
    envelopes = Vector{ForceEnvelopes}()
    forceresults = Vector{Vector{InternalForces}}()

    for load in loads
        solve!(model, load)
        push!(forceresults, forces(model, increment))
    end

    n_elem = length(first(forceresults))
    n_loads = length(forceresults)
    sizehint!(envelopes, n_elem)

    for i = 1:n_elem
        ref = first(forceresults)[i]
        e   = ref.element
        res = ref.resolution
        x   = ref.x
        n_pts = length(x)

        # Preallocate matrices (pts × load cases) — avoids hcat+getproperty broadcasts
        Pm  = Matrix{Float64}(undef, n_pts, n_loads)
        Mym = Matrix{Float64}(undef, n_pts, n_loads)
        Vym = Matrix{Float64}(undef, n_pts, n_loads)
        Mzm = Matrix{Float64}(undef, n_pts, n_loads)
        Vzm = Matrix{Float64}(undef, n_pts, n_loads)

        @inbounds for j in 1:n_loads
            fr = forceresults[j][i]
            for k in 1:n_pts
                Pm[k, j]  = fr.P[k]
                Mym[k, j] = fr.My[k]
                Vym[k, j] = fr.Vy[k]
                Mzm[k, j] = fr.Mz[k]
                Vzm[k, j] = fr.Vz[k]
            end
        end

        # Compute envelopes row-by-row
        Plow  = Vector{Float64}(undef, n_pts);  Phigh  = Vector{Float64}(undef, n_pts)
        Mylow = Vector{Float64}(undef, n_pts);  Myhigh = Vector{Float64}(undef, n_pts)
        Vylow = Vector{Float64}(undef, n_pts);  Vyhigh = Vector{Float64}(undef, n_pts)
        Mzlow = Vector{Float64}(undef, n_pts);  Mzhigh = Vector{Float64}(undef, n_pts)
        Vzlow = Vector{Float64}(undef, n_pts);  Vzhigh = Vector{Float64}(undef, n_pts)

        @inbounds for k in 1:n_pts
            lo, hi = Pm[k, 1], Pm[k, 1]
            for j in 2:n_loads; v = Pm[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Plow[k] = lo; Phigh[k] = hi

            lo, hi = Mym[k, 1], Mym[k, 1]
            for j in 2:n_loads; v = Mym[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Mylow[k] = lo; Myhigh[k] = hi

            lo, hi = Vym[k, 1], Vym[k, 1]
            for j in 2:n_loads; v = Vym[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Vylow[k] = lo; Vyhigh[k] = hi

            lo, hi = Mzm[k, 1], Mzm[k, 1]
            for j in 2:n_loads; v = Mzm[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Mzlow[k] = lo; Mzhigh[k] = hi

            lo, hi = Vzm[k, 1], Vzm[k, 1]
            for j in 2:n_loads; v = Vzm[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Vzlow[k] = lo; Vzhigh[k] = hi
        end

        push!(envelopes, ForceEnvelopes(e, res, x,
            Plow, Phigh, Mylow, Myhigh, Vylow, Vyhigh, Mzlow, Mzhigh, Vzlow, Vzhigh))
    end

    return envelopes
end

# =============================================================================
# Convenience: moments! / shear! for frame elements
# =============================================================================

"""
    moments!(My, Mz, element, model; resolution=20)

Write bending moments My (about local Y) and Mz (about local Z) into
pre-allocated vectors.  Returns `(My, Mz)`.
"""
function moments!(My::Vector{Float64}, Mz::Vector{Float64},
                  element::AbstractElement, model::Union{FrameModel, Model};
                  resolution::Int = 20)
    eif = ElementInternalForces(element, model; resolution=resolution)
    copyto!(My, eif.My)
    copyto!(Mz, eif.Mz)
    return (My, Mz)
end

"""
    shears!(Vy, Vz, element, model; resolution=20)

Write shear forces Vy and Vz into pre-allocated vectors.  Returns `(Vy, Vz)`.
"""
function shears!(Vy::Vector{Float64}, Vz::Vector{Float64},
                 element::AbstractElement, model::Union{FrameModel, Model};
                 resolution::Int = 20)
    eif = ElementInternalForces(element, model; resolution=resolution)
    copyto!(Vy, eif.Vy)
    copyto!(Vz, eif.Vz)
    return (Vy, Vz)
end
