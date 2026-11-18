# =============================================================================
# accumulate_force!  —  TributaryLoad
# =============================================================================

"""
    _tributary_end_moments(element, mA, mB) -> (eA, eB)

End moments of a load whose fixed-fixed end moment magnitudes are `mA` and
`mB`, after the element's releases. The carry-over is the one
[`q`](@ref) applies to the fixed-end forces, so the recovered field and the
assembled load vector describe the same member.
"""
_tributary_end_moments(::Element{FixedFixed}, mA, mB) = (mA, mB)
_tributary_end_moments(::Element{FreeFixed}, mA, mB) = (0.0, mB + mA / 2)
_tributary_end_moments(::Element{FixedFree}, mA, mB) = (mA + mB / 2, 0.0)
_tributary_end_moments(::Element{FreeFree}, mA, mB) = (0.0, 0.0)
_tributary_end_moments(::Element{Joist}, mA, mB) = (0.0, 0.0)

"""
    _segment_partial_integrals(a, b, wa, wb, x) -> (F, Fm)

For `w(ξ)` linear from `wa` at `a` to `wb` at `b`, the integrals up to station
`x`: `F = ∫ w dξ` and `Fm = ∫ w (x - ξ) dξ`, both over `[a, min(x, b)]`.
"""
@inline function _segment_partial_integrals(a, b, wa, wb, x)
    x <= a && return (0.0, 0.0)
    k = (wb - wa) / (b - a)
    T = min(x, b) - a
    D = x - a
    F = wa * T + k * T^2 / 2
    Fm = wa * (D * T - T^2 / 2) + k * (D * T^2 / 2 - T^3 / 3)
    return (F, Fm)
end

"""
Accumulate the internal forces caused by a TributaryLoad to the current element.

The load is piecewise linear in intensity. Each segment is integrated exactly:
the simple-span field `Ra x - ∫ w (x - ξ) dξ`, minus the end moments of the
element's release, interpolated linearly. The fixed-fixed end moments are the
same integrals as the fixed-end forces ([`_fem_uniform_partial`](@ref),
[`_fem_triangular_partial`](@ref)). Sign conventions match the `LineLoad`
method: a uniform tributary load gives the same field as the equal `LineLoad`.
"""
function accumulate_force!(load::TributaryLoad,
    xvals::AbstractVector{Float64},
    P::Vector{Float64},
    My::Vector{Float64},
    Vy::Vector{Float64},
    Mz::Vector{Float64},
    Vz::Vector{Float64},
    ep::ElementProps)

    L = ep.L
    d1, d2, d3 = load.direction
    ux, uy, uz = _rotate_to_local(ep.Rv, d1, d2, d3)
    w = intensities(load)
    positions = load.positions

    # Start reaction of the simple span and fixed-fixed end moment magnitudes,
    # per unit load direction.
    Ra = 0.0
    mA = 0.0
    mB = 0.0
    for i in 1:(length(positions) - 1)
        a = positions[i] * L
        b = positions[i + 1] * L
        ℓ = b - a
        ℓ < 1e-12 && continue
        wa, wb = w[i], w[i + 1]
        k = (wb - wa) / ℓ
        Ra += (wa * ℓ + k * ℓ^2 / 2) -
              (wa * (a * ℓ + ℓ^2 / 2) + k * (a * ℓ^2 / 2 + ℓ^3 / 3)) / L
        u = _fem_uniform_partial(L, a, b, wa)
        t = _fem_triangular_partial(L, a, b, wb - wa)
        mA -= u.M_A + t.M_A
        mB += u.M_B + t.M_B
    end
    eA, eB = _tributary_end_moments(load.element, mA, mB)
    Vend = (eA - eB) / L

    @inbounds for (j, x) in enumerate(xvals)
        F = 0.0
        Fm = 0.0
        for i in 1:(length(positions) - 1)
            a = positions[i] * L
            b = positions[i + 1] * L
            (b - a) < 1e-12 && continue
            Fi, Fmi = _segment_partial_integrals(a, b, w[i], w[i + 1], x)
            F += Fi
            Fm += Fmi
        end
        M = Ra * x - Fm - (eA * (1 - x / L) + eB * x / L)
        V = Ra - F + Vend
        P[j]  += ux * (Ra - F)
        My[j] += uy * M
        Vy[j] += uy * V
        Mz[j] += uz * M
        Vz[j] += uz * V
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
