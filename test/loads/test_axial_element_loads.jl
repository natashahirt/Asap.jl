using Asap
using Test
using Unitful

# Axial components of element loads on a bar of uniform EA fixed axially at
# both ends. Statics + compatibility: a point load P at distance a from the
# start sends P b / L to the start and P a / L to the end (b = L - a); a
# distributed load w(x) sends ∫ w (L - x) / L and ∫ w x / L. Internal axial
# force is tension positive.

const _AX_SEC = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")

_ax_node(x, fix) = Node([x * u"m", 0.0u"m", 0.0u"m"], fix)

"""Axial force at station `x` (m) of a solved element, by internal-force recovery."""
function _axial_at(el, model, x)
    f = ElementInternalForces(el, model; resolution = 401)
    i = argmin(abs.(f.x .- x))
    @assert abs(f.x[i] - x) < 1e-9
    return f.P[i]
end

@testset "Axial element loads" begin
    @testset "PointLoad axial split follows position" begin
        # Two collinear elements, outer ends fixed, middle node free. P at
        # quarter of element 1 = L/8 globally, so the exact field is
        # +P(7/8) on [0, L/8] and -P/8 on [L/8, L]. The old P/2 split gave
        # -P/4 in element 2.
        L = 8.0
        P = 100e3
        n1 = _ax_node(0.0, :fixed); n2 = _ax_node(L / 2, :free); n3 = _ax_node(L, :fixed)
        e1 = Element(n1, n2, _AX_SEC); e2 = Element(n2, n3, _AX_SEC)
        ld = PointLoad(e1, 0.25, [P, 0.0, 0.0] .* u"N")
        model = Model([n1, n2, n3], [e1, e2], AbstractLoad[ld])
        solve!(model)

        @test _axial_at(e2, model, 0.0) ≈ -P / 8 rtol = 1e-9
        @test _axial_at(e2, model, L / 2) ≈ -P / 8 rtol = 1e-9
        @test _axial_at(e1, model, 0.0) ≈ 7P / 8 rtol = 1e-9
        @test _axial_at(e1, model, 0.5) ≈ 7P / 8 rtol = 1e-9     # x < a = 1 m
        @test _axial_at(e1, model, 2.0) ≈ -P / 8 rtol = 1e-9     # x > a
        @test _axial_at(e1, model, L / 2) ≈ -P / 8 rtol = 1e-9
        # Support reactions balance the load with the same split.
        @test ustrip(u"N", n1.reaction[1]) ≈ -7P / 8 rtol = 1e-9
        @test ustrip(u"N", n3.reaction[1]) ≈ -P / 8 rtol = 1e-9
    end

    @testset "PointLoad axial on a single fixed-fixed element" begin
        # Ku = 0, so the recovered field is the load term alone: P b/L, -P a/L.
        L = 4.0; P = 100e3
        n1 = _ax_node(0.0, :fixed); n2 = _ax_node(L, :fixed)
        e = Element(n1, n2, _AX_SEC)
        model = Model([n1, n2], [e], AbstractLoad[PointLoad(e, 0.25, [P, 0.0, 0.0] .* u"N")])
        solve!(model)
        @test _axial_at(e, model, 0.5) ≈ 0.75P rtol = 1e-9
        @test _axial_at(e, model, 3.0) ≈ -0.25P rtol = 1e-9
        @test ustrip(u"N", n1.reaction[1]) ≈ -0.75P rtol = 1e-9
        @test ustrip(u"N", n2.reaction[1]) ≈ -0.25P rtol = 1e-9
    end

    @testset "Triangular axial TributaryLoad" begin
        # Single element, L = 4.5 m, w = 0 → 20 kN/m along +x.
        # N_A = ∫ w (L-x)/L = wmax L / 6 = 15 kN, N_B = wmax L / 3 = 30 kN.
        # N(L) = N_A - wmax L / 2 = -30 kN. The equal split gave -22.5 kN.
        L = 4.5; wmax = 20e3
        n1 = _ax_node(0.0, :fixed); n2 = _ax_node(L, :fixed)
        e = Element(n1, n2, _AX_SEC)
        ld = TributaryLoad(e, [0.0, 1.0], [0.0u"m", 1.0u"m"], wmax * u"Pa", (1.0, 0.0, 0.0))
        model = Model([n1, n2], [e], AbstractLoad[ld])
        solve!(model)
        @test _axial_at(e, model, 0.0) ≈ 15e3 rtol = 1e-9
        @test _axial_at(e, model, L) ≈ -30e3 rtol = 1e-9
        @test ustrip(u"N", n1.reaction[1]) ≈ -15e3 rtol = 1e-9
        @test ustrip(u"N", n2.reaction[1]) ≈ -30e3 rtol = 1e-9

        # Same load on element 1 of two collinear elements (fixed ends, free
        # middle): the fixed-end forces now drive the solution. Global L = 9,
        # load on [0, 4.5]: N_end = ∫ w x / L = wmax (L/2)^2 / 3 / L = 15 kN,
        # so element 2 carries -15 kN and the start carries +30 kN.
        L = 9.0
        n1 = _ax_node(0.0, :fixed); n2 = _ax_node(L / 2, :free); n3 = _ax_node(L, :fixed)
        e1 = Element(n1, n2, _AX_SEC); e2 = Element(n2, n3, _AX_SEC)
        ld = TributaryLoad(e1, [0.0, 1.0], [0.0u"m", 1.0u"m"], wmax * u"Pa", (1.0, 0.0, 0.0))
        model = Model([n1, n2, n3], [e1, e2], AbstractLoad[ld])
        solve!(model)
        @test _axial_at(e2, model, 0.0) ≈ -15e3 rtol = 1e-9
        @test _axial_at(e1, model, 0.0) ≈ 30e3 rtol = 1e-9
        @test _axial_at(e1, model, L / 2) ≈ -15e3 rtol = 1e-9
    end

    @testset "_fem_axial_partial closed forms" begin
        # Uniform w on [a, b]: N_B = w (b² - a²) / (2L).
        NA, NB = Asap._fem_axial_partial(8.0, 1.0, 5.0, 250.0, 0.0, 0.0)
        @test NB ≈ 250.0 * (25.0 - 1.0) / 16.0 rtol = 1e-12
        @test NA + NB ≈ 1000.0 rtol = 1e-12
        # w = τ² on [0, L]: total L/3, N_B = ∫ x³/L³ dx = L/4.
        NA, NB = Asap._fem_axial_partial(2.0, 0.0, 2.0, 0.0, 0.0, 1.0)
        @test NB ≈ 0.5 rtol = 1e-12
        @test NA ≈ 2.0 / 3 - 0.5 rtol = 1e-12
        # The quadratic and linear accumulators agree on a linear axial load.
        fq = zeros(12); fl = zeros(12)
        Asap._fem_quadratic_segment!(fq, 8.0, 0.25, 0.75, 100.0, 50.0, 0.0, -1.0, 0.0, 0.0)
        Asap._fem_linear_segment!(fl, 8.0, 0.25, 0.75, 100.0, 150.0, -1.0, 0.0, 0.0)
        @test fq ≈ fl rtol = 1e-12
        # 100 → 150 N/m on [2, 6] of L = 8: total 500 N,
        # ∫ w x = ℓ/6 (w_a (2a + b) + w_b (a + 2b)) = 6200/3, so N_B = 775/3.
        @test fl[7] ≈ -775 / 3 rtol = 1e-12
        @test fl[1] ≈ -(500 - 775 / 3) rtol = 1e-12
    end

    @testset "ThetaLine_fixedfixed is the slope of DLine_fixedfixed" begin
        w = 1.5e3; L = 6.0; E = 200e9; I = 3e-5
        h = 1e-5
        for x in (0.3, 1.0, 2.2, 3.7, 5.1)
            fd = (Asap.DLine_fixedfixed(w, L, x + h, E, I) -
                  Asap.DLine_fixedfixed(w, L, x - h, E, I)) / 2h
            @test Asap.ThetaLine_fixedfixed(w, L, x, E, I) ≈ fd rtol = 1e-7
        end
        @test Asap.ThetaLine_fixedfixed(w, L, L / 2, E, I) == 0.0
        @test Asap.ThetaLine_fixedfixed(w, L, 0.0, E, I) == 0.0
        @test Asap.ThetaLine_fixedfixed(w, L, L, E, I) == 0.0
        # Closed form at L/4: w x (L-x)(L-2x) / (12 E I).
        @test Asap.ThetaLine_fixedfixed(w, L, L / 4, E, I) ≈
              w * (L / 4) * (3L / 4) * (L / 2) / (12E * I) rtol = 1e-14
    end
end
