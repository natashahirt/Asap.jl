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


