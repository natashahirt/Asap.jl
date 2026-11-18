using Asap
using Test
using Unitful

# Vertical service deflection from a displacement vector and the case loads
# (`element_deflection_profile` / `element_max_deflections`).
#
# References (standard beam tables, uniform load w on span L):
#   - simple span:     δ_mid = 5 w L⁴ / (384 E I)
#   - fixed-fixed:     δ_mid =   w L⁴ / (384 E I)
# With the beam rolled so that neither principal axis is vertical, a vertical
# load splits into w·(ly·ẑ) along local y and w·(lz·ẑ) along local z. The two
# planes bend independently (I about local z is Ix, about local y is Iy), so the
# vertical deflection is the table value with 1/I replaced by
#   (ly·ẑ)²/Ix + (lz·ẑ)²/Iy.
#
# Each span is ONE element, as production models one element per skeleton
# edge. Before the fix the nodal-only interpolation of local z returned 0 for
# a simple span (the 9 m, 20 kN/m case below reported ~0 instead of 40.24 mm).

const _VD_E = 200e9
const _VD_IX = 8.0e-5
const _VD_IY = 2.0e-5
const _VD_SEC = Section(0.01u"m^2", _VD_E * u"Pa", 77e9u"Pa",
                        _VD_IX * u"m^4", _VD_IY * u"m^4", 1e-5u"m^4")

"""Fixity that pins translations and the torsion rotation, freeing bending."""
function _vd_pin(axis)
    dofs = [false, false, false, true, true, true]
    dofs[abs(axis[1]) > 0.5 ? 4 : 5] = false
    return dofs
end

"""One-element beam from the origin along `dir`; returns (model, element, loads)."""
function _vd_beam(L, w; dir = (1.0, 0.0, 0.0), release = :freefree, supports = :pinned,
                  load = :tributary, roll = nothing)
    p1 = [0.0u"m", 0.0u"m", 3.0u"m"]
    p2 = [dir[1] * L * u"m", dir[2] * L * u"m", 3.0u"m"]
    fix = if supports === :fixed
        (:fixed, :fixed)
    elseif release === :freefree || release === :joist
        # A released element puts no bending stiffness on the node, so the
        # node rotations are restrained; the element still rotates freely.
        (:fixed, :fixed)
    else
        (_vd_pin(dir), _vd_pin(dir))
    end
    n1 = Node(p1, fix[1])
    n2 = Node(p2, fix[2])
    el = Element(n1, n2, _VD_SEC; release = release)
    roll === nothing || (el.Ψ = roll)
    ld = load === :line ?
        LineLoad(el, [0.0u"N/m", 0.0u"N/m", -w * u"N/m"]) :
        TributaryLoad(el, [0.0, 1.0], [1.0u"m", 1.0u"m"], w * u"Pa", (0.0, 0.0, -1.0))
    loads = Asap.AbstractLoad[ld]
    model = Model([n1, n2], Asap.FrameElement[el], loads)
    solve!(model)
    return model, el, loads
end

"""Vertical-plane flexibility (ly·ẑ)²/Ix + (lz·ẑ)²/Iy of a horizontal element."""
_vd_flex(el) = el.LCS[2][3]^2 / _VD_IX + el.LCS[3][3]^2 / _VD_IY

@testset "vertical deflection of one-element spans" begin
    L = 9.0
    w = 20_000.0

    @testset "simple span, $release element, $load load, along $dir" for
            release in (:freefree, :joist, :fixedfixed),
            load in (:line, :tributary),
            dir in ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0))
        model, el, loads = _vd_beam(L, w; dir, release, load)
        δ_exact = 5 * w * L^4 / 384 * _vd_flex(el) / _VD_E
        prof = Asap.element_deflection_profile(el, model.u, loads; resolution = 21)
        # Station 11 of 21 is midspan; downward is negative.
        @test prof[11] ≈ -δ_exact rtol = 1e-6
        @test abs(prof[1]) < 1e-12 && abs(prof[end]) < 1e-12
        δmax = Asap.element_max_deflections(model, model.u, loads; resolution = 21)
        @test δmax[el.elementID] ≈ δ_exact rtol = 1e-6
    end

    @testset "audit case: 9 m, 20 kN/m, default roll reads the strong axis" begin
        model, el, loads = _vd_beam(L, w)
        @test Asap.element_max_deflections(model, model.u, loads; resolution = 21)[el.elementID] ≈
              5 * w * L^4 / (384 * _VD_E * _VD_IX) rtol = 1e-6
        @test 5 * w * L^4 / (384 * _VD_E * _VD_IX) ≈ 0.10679 rtol = 1e-4
    end

    @testset "fixed-fixed supports, $load load, along $dir" for
            load in (:line, :tributary),
            dir in ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0))
        model, el, loads = _vd_beam(L, w; dir, release = :fixedfixed, supports = :fixed, load)
        δ_exact = w * L^4 / 384 * _vd_flex(el) / _VD_E
        prof = Asap.element_deflection_profile(el, model.u, loads; resolution = 21)
        @test prof[11] ≈ -δ_exact rtol = 1e-6
    end

    @testset "roll Ψ = $roll" for roll in (0.0, 0.7, pi / 2)
        model, el, loads = _vd_beam(L, w; release = :fixedfixed, roll = roll)
        # Ψ = 0 puts local z vertical (weak axis Iy); π/2 puts local y vertical.
        roll == 0.0 && @test abs(el.LCS[3][3]) ≈ 1.0
        roll == pi / 2 && @test abs(el.LCS[2][3]) ≈ 1.0
        δ_exact = 5 * w * L^4 / 384 * _vd_flex(el) / _VD_E
        prof = Asap.element_deflection_profile(el, model.u, loads; resolution = 21)
        @test prof[11] ≈ -δ_exact rtol = 1e-6
    end

    @testset "loads of another element are refused" begin
        model, el, loads = _vd_beam(L, w)
        _, other, other_loads = _vd_beam(L, w)
        @test_throws ArgumentError Asap.element_deflection_profile(el, model.u, other_loads)
        @test_throws DimensionMismatch Asap.element_max_deflections(model, zeros(3), loads)
    end
end
