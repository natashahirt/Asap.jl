using Asap
using Test
using Unitful

# TributaryLoad on every end release: fixed-end forces, internal forces, and
# deflections.
#
# References:
#   - a uniform TributaryLoad against the equal LineLoad (closed-form MLine,
#     VLine, PLine, DLine for each release);
#   - a varying TributaryLoad on fixed supports against closed-form fields:
#     M(x) = Ra x - ∫ w (x - ξ) dξ - [eA (1 - x/L) + eB x/L], with
#     eA = ∫ w x (L - x)² / L², eB = ∫ w x² (L - x) / L² for fixed-fixed and
#     the carry-over eB + eA/2 (or eA + eB/2) for one released end. Triangle,
#     0 → w: fixed-fixed eA = wL²/30, eB = wL²/20; pinned-fixed eB = wL²/15;
#     fixed-pinned eA = 7wL²/120 (standard beam tables);
#   - a varying TributaryLoad on flexible columns against the same beam cut
#     into 64 elements, each with a uniform LineLoad at the local mean.

const _TR_SEC = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 1e-4u"m^4", 2e-4u"m^4", 1e-4u"m^4")
const _TR_COL = Section(0.02u"m^2", 200e9u"Pa", 77e9u"Pa", 5e-5u"m^4", 5e-5u"m^4", 5e-5u"m^4")
const _TR_W = 1000.0
const _TR_L = 10.0
const _TR_RELEASES = (:fixedfixed, :freefixed, :fixedfree, :freefree, :joist)

_tr_el(n1, n2, rel) = Element(n1, n2, _TR_SEC; release = rel)

function _tr_fields(model, el; resolution = 9)
    loads = get(Asap.get_elemental_loads(model), el.elementID, Asap.AbstractLoad[])
    r = Asap.ElementForceAndDisplacement(el, loads; resolution = resolution)
    return (P = r.forces.P, My = r.forces.My, Vy = r.forces.Vy,
            Mz = r.forces.Mz, Vz = r.forces.Vz, D = r.displacements.uglobal)
end

_tr_rel_err(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), 1e-9)

@testset "TributaryLoad fields on every release" begin
    @testset "uniform load equals LineLoad ($rel, $where)" for rel in _TR_RELEASES,
            where in (:fixed_supports, :columns, :sloped)
        function build(kind)
            z1 = where === :sloped ? 2.0 : (where === :columns ? 3.0 : 0.0)
            z0 = where === :columns ? 3.0 : 0.0
            fix = where === :fixed_supports ? :fixed : :free
            n1 = Node([0.0u"m", 0.0u"m", z0 * u"m"], fix)
            n2 = Node([_TR_L * u"m", 0.0u"m", z1 * u"m"], fix)
            nodes = [n1, n2]
            els = Asap.FrameElement[]
            el = _tr_el(n1, n2, rel)
            push!(els, el)
            if where !== :fixed_supports
                g1 = Node([0.0u"m", 0.0u"m", (z0 - 3.0) * u"m"], :fixed)
                g2 = Node([_TR_L * u"m", 0.0u"m", (z1 - 3.0) * u"m"], :fixed)
                append!(nodes, [g1, g2])
                append!(els, [Element(g1, n1, _TR_COL), Element(g2, n2, _TR_COL)])
            end
            Lel = ustrip(u"m", el.length)
            # The same total load per unit element length in both models.
            ld = kind === :line ?
                LineLoad(el, [0.0u"N/m", 0.0u"N/m", -_TR_W * u"N/m"]) :
                TributaryLoad(el, [0.0, 1.0], [1.0u"m", 1.0u"m"], _TR_W * u"Pa", (0.0, 0.0, -1.0))
            m = Model(nodes, els, [ld])
            solve!(m)
            return (m = m, el = el, L = Lel)
        end
        line = build(:line)
        trib = build(:trib)
        @test trib.m.reactions ≈ line.m.reactions atol = 1e-6 * _TR_W * _TR_L^2
        fl = _tr_fields(line.m, line.el)
        ft = _tr_fields(trib.m, trib.el)
        for k in (:P, :My, :Vy, :Mz, :Vz, :D)
            @test getproperty(ft, k) ≈ getproperty(fl, k) atol = 1e-8 * max(1.0, maximum(abs.(getproperty(fl, k))))
        end
    end

    @testset "varying load on fixed supports ($rel, $shape)" for rel in _TR_RELEASES,
            (shape, c0, c1) in (("triangle", 0.0, 1.0), ("trapezoid", 0.5, 1.5))
        L = _TR_L
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L * u"m", 0.0u"m", 0.0u"m"], :fixed)
        el = _tr_el(n1, n2, rel)
        ld = TributaryLoad(el, [0.0, 1.0], [c0 * u"m", (c0 + c1) * u"m"], _TR_W * u"Pa", (0.0, 0.0, -1.0))
        m = Model([n1, n2], [el], [ld]); solve!(m)
        f = _tr_fields(m, el)
        # Closed form for w(x) = W (c0 + c1 x/L).
        W = _TR_W
        tot = W * (c0 + c1 / 2) * L
        Ra = W * L * (c0 / 2 + c1 / 6)
        mA = W * L^2 * (c0 / 12 + c1 / 30)
        mB = W * L^2 * (c0 / 12 + c1 / 20)
        eA, eB = rel === :fixedfixed ? (mA, mB) : rel === :freefixed ? (0.0, mB + mA / 2) :
                 rel === :fixedfree ? (mA + mB / 2, 0.0) : (0.0, 0.0)
        xs = range(0, L, 9)
        Mss(x) = Ra * x - W * (c0 * x^2 / 2 + c1 * x^3 / (6L))
        M_hand = [Mss(x) - eA * (1 - x / L) - eB * x / L for x in xs]
        @test f.My ≈ M_hand atol = 1e-8 * W * L^2
        V_hand = [Ra - W * (c0 * x + c1 * x^2 / (2L)) + (eA - eB) / L for x in xs]
        @test f.Vy ≈ V_hand atol = 1e-8 * W * L
        @test sum(m.reactions[[3, 9]]) ≈ tot rtol = 1e-10
        if shape == "triangle"
            # The hand end moments against the beam-table values.
            if rel === :fixedfixed
                @test eA ≈ W * L^2 / 30
                @test eB ≈ W * L^2 / 20
            elseif rel === :freefixed
                @test eB ≈ W * L^2 / 15
            elseif rel === :fixedfree
                @test eA ≈ 7 * W * L^2 / 120
            end
            # End moments read off the FE field.
            @test -f.My[1] ≈ eA atol = 1e-8 * W * L^2
            @test -f.My[end] ≈ eB atol = 1e-8 * W * L^2
        end
    end

    # A beam released at both ends on columns has no torsional restraint; the
    # fixed-support cases above cover those releases.
    @testset "varying load on columns matches a 64-element LineLoad beam ($rel, $shape)" for rel in (:fixedfixed, :freefixed, :fixedfree),
            (shape, positions, widths) in (("triangle", [0.0, 1.0], [0.0, 1.0]),
                                           ("four-piece", [0.0, 0.2, 0.55, 1.0], [0.3, 1.2, 0.8, 0.0]))
        L = _TR_L
        wfun(s) = _TR_W * (let i = clamp(searchsortedlast(positions, s), 1, length(positions) - 1)
            widths[i] + (widths[i + 1] - widths[i]) * (s - positions[i]) / (positions[i + 1] - positions[i])
        end)
        function frame(n, trib)
            top(x) = Node([x * u"m", 0.0u"m", 3.0u"m"], :free)
            chain = [top(L * k / n) for k in 0:n]
            gA = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            gB = Node([L * u"m", 0.0u"m", 0.0u"m"], :fixed)
            els = Asap.FrameElement[]; lds = Asap.AbstractLoad[]
            start_rel = rel in (:freefixed, :freefree, :joist)
            end_rel = rel in (:fixedfree, :freefree, :joist)
            for k in 1:n
                r = n == 1 ? rel :
                    k == 1 && start_rel ? :freefixed :
                    k == n && end_rel ? :fixedfree : :fixedfixed
                e = _tr_el(chain[k], chain[k + 1], r); push!(els, e)
                if trib
                    push!(lds, TributaryLoad(e, positions, widths .* u"m", _TR_W * u"Pa", (0.0, 0.0, -1.0)))
                else
                    s0 = (k - 1) / n; s1 = k / n
                    wm = (wfun(s0) + 4wfun((s0 + s1) / 2) + wfun(s1)) / 6
                    push!(lds, LineLoad(e, [0.0u"N/m", 0.0u"N/m", -wm * u"N/m"]))
                end
            end
            append!(els, [Element(gA, chain[1], _TR_COL), Element(gB, chain[end], _TR_COL)])
            m = Model([chain; gA; gB], els, lds); solve!(m)
            return (m = m, els = els, chain = chain)
        end
        t = frame(1, true)
        r = frame(64, false)
        ft = _tr_fields(t.m, t.els[1]; resolution = 9)
        stations = 0:8:64
        My_ref = Float64[]
        for j in stations
            # Station j is the start of element j + 1, or the end of element 64.
            k = j == 64 ? 64 : j + 1
            push!(My_ref, _tr_fields(r.m, r.els[k]; resolution = 2).My[j == 64 ? 2 : 1])
        end
        @test _tr_rel_err(ft.My, My_ref) < 2e-3
        uz_ref = [ustrip(u"m", r.chain[j + 1].displacement[3]) for j in stations]
        @test _tr_rel_err(ft.D[3, :], uz_ref) < 2e-3
    end

    @testset "release_end_rotations! gives zero curvature at a released end" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
        L = 4.0
        v1, v2, θ1, θ2 = 0.002, -0.003, 0.0011, -0.0007
        # Hermite curvature at ξ = 0 and ξ = 1 for (v1, s1, v2, s2).
        k0(v1, s1, v2, s2) = 6(v2 - v1) / L^2 - (4s1 + 2s2) / L
        k1(v1, s1, v2, s2) = -6(v2 - v1) / L^2 + (2s1 + 4s2) / L
        for (rel, zero_start, zero_end) in ((:freefixed, true, false), (:fixedfree, false, true),
                                            (:freefree, true, true), (:joist, true, true))
            el = Element(n1, n2, _TR_SEC; release = rel)
            ul = zeros(12)
            ul[2], ul[6], ul[8], ul[12] = v1, θ1, v2, θ2
            ul[3], ul[5], ul[9], ul[11] = v1, -θ1, v2, -θ2
            ul .*= etype2DOF[typeof(el)]
            Asap.release_end_rotations!(ul, el, L)
            sy = (ul[6], ul[12]); sz = (-ul[5], -ul[11])
            for s in (sy, sz)
                zero_start && @test k0(v1, s[1], v2, s[2]) ≈ 0 atol = 1e-12
                zero_end && @test k1(v1, s[1], v2, s[2]) ≈ 0 atol = 1e-12
                !zero_start && @test s[1] == θ1
                !zero_end && @test s[2] == θ2
            end
        end
    end
end
