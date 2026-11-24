using Asap
using Test
using Unitful
using LinearAlgebra: dot, norm

# Virtual-work drift participation: δ at a free DOF decomposes into per-element
# coefficients c_e = u_vᵀ K_e u_r, with Σ c_e == δ (Maxwell–Betti). These tests
# check (1) the exact decomposition identity, (2) absolute magnitudes against
# hand calculations (axial PL/EA, cantilever PL³/3EI), (3) a portal-frame sway
# where the columns dominate participation, and (4) that the virtual solve
# reuses the cached factorization (no extra factorization).

@testset "drift participation" begin

    @testset "axial bar: single-element participation equals PL/EA" begin
        L = 2.0
        A = 1.0e-2
        E = 200.0e9
        P = 1000.0
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L * 1.0u"m", 0.0u"m", 0.0u"m"], :xfree)
        sec = Section(A * u"m^2", E * u"Pa", 80.0e9u"Pa",
                      1.0e-5u"m^4", 1.0e-5u"m^4", 1.0e-6u"m^4")
        el = Element(n1, n2, sec)
        model = FrameModel([n1, n2], [el],
                           [NodeForce(n2, [P * 1.0u"N", 0.0u"N", 0.0u"N"])])
        solve!(model)

        dof = n2.globalID[1]
        dp = drift_participation(model, dof)

        δ_exact = P * L / (E * A)                 # PL/EA
        @test dp.δ ≈ δ_exact rtol = 1e-9
        @test length(dp.c) == 1
        @test dp.c[1] ≈ δ_exact rtol = 1e-9       # the one element carries all of δ
        @test abs(dp.nonframe) < 1e-12 * max(1.0, abs(dp.δ))
        # Node + direction overload resolves to the same DOF.
        @test drift_participation(model, n2, 1).δ ≈ dp.δ
    end

    @testset "cantilever: Σ c == δ and δ ≈ PL³/3EI" begin
        # Slender cantilever along X, transverse tip load in Z, split into
        # several elements so the decomposition is non-trivial. Equal Ix = Iy
        # makes the bending inertia unambiguous; slenderness makes shear < 1%.
        Ltot = 10.0
        nseg = 5
        I = 1.0e-5
        A = 1.0e-2
        E = 200.0e9
        P = 500.0
        sec = Section(A * u"m^2", E * u"Pa", 80.0e9u"Pa",
                      I * u"m^4", I * u"m^4", 1.0e-6u"m^4")
        nodes = [Node([x * 1.0u"m", 0.0u"m", 0.0u"m"],
                      x == 0.0 ? :fixed : :free)
                 for x in range(0.0, Ltot; length = nseg + 1)]
        els = [Element(nodes[i], nodes[i + 1], sec) for i in 1:nseg]
        tip = nodes[end]
        model = FrameModel(nodes, els,
                           [NodeForce(tip, [0.0u"N", 0.0u"N", -P * 1.0u"N"])])
        solve!(model)

        dof = tip.globalID[3]                      # local Z translation at the tip
        dp = drift_participation(model, dof)

        # (1) Exact decomposition: Σ c_e == δ == model.u[dof], no non-frame share.
        @test sum(dp.c) ≈ dp.δ rtol = 1e-9
        @test dp.δ ≈ model.u[dof] rtol = 1e-12
        @test abs(dp.nonframe) < 1e-8 * max(1.0, abs(dp.δ))
        @test length(dp.c) == nseg

        # (2) Magnitude against Euler–Bernoulli PL³/3EI (slender ⇒ shear < ~1%).
        δ_beam = -P * Ltot^3 / (3 * E * I)         # negative: load is −Z
        @test dp.δ ≈ δ_beam rtol = 2e-2
        # Every segment participates with the same sign as the drift.
        @test all(c -> c * dp.δ ≥ -1e-9, dp.c)
    end

    @testset "portal frame sway: columns dominate participation" begin
        H = 3.0
        W = 4.0
        E = 200.0e9
        P = 10_000.0
        col = Section(1.0e-2u"m^2", E * u"Pa", 80.0e9u"Pa",
                      8.0e-5u"m^4", 8.0e-5u"m^4", 1.0e-6u"m^4")
        beam = Section(1.2e-2u"m^2", E * u"Pa", 80.0e9u"Pa",
                       1.5e-4u"m^4", 1.5e-4u"m^4", 1.0e-6u"m^4")
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)     # base, column 1
        n2 = Node([0.0u"m", 0.0u"m", H * 1.0u"m"], :free)  # top, column 1
        n3 = Node([W * 1.0u"m", 0.0u"m", H * 1.0u"m"], :free)
        n4 = Node([W * 1.0u"m", 0.0u"m", 0.0u"m"], :fixed) # base, column 2
        c1 = Element(n1, n2, col)                          # columns bend in X–Z
        c2 = Element(n4, n3, col)
        gb = Element(n2, n3, beam)                         # girder
        model = FrameModel([n1, n2, n3, n4], [c1, c2, gb],
                           [NodeForce(n2, [P * 1.0u"N", 0.0u"N", 0.0u"N"])])
        solve!(model)

        dof = n2.globalID[1]                                # X sway at the loaded top
        dp = drift_participation(model, dof)

        @test dp.δ > 0
        @test sum(dp.c) ≈ dp.δ rtol = 1e-9
        @test abs(dp.nonframe) < 1e-8 * max(1.0, abs(dp.δ))

        # Columns (ids 1,2) carry far more drift participation than the girder (3).
        cbyid = Dict(zip(dp.element_ids, dp.c))
        col_share = cbyid[c1.elementID] + cbyid[c2.elementID]
        @test col_share > 0.5 * dp.δ
    end

    @testset "virtual solve reuses the cached factorization" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([3.0u"m", 0.0u"m", 0.0u"m"], :free)
        sec = Section(1.0e-2u"m^2", 200.0e9u"Pa", 80.0e9u"Pa",
                      5.0e-5u"m^4", 5.0e-5u"m^4", 1.0e-6u"m^4")
        el = Element(n1, n2, sec)
        model = FrameModel([n1, n2], [el],
                           [NodeForce(n2, [0.0u"N", 0.0u"N", -1000.0u"N"])])
        solve!(model)

        before = Asap._N_FACTORIZATIONS[]
        dp = drift_participation(model, n2.globalID[3])
        dp2 = drift_participation(model, n2.globalID[1])   # a second virtual solve
        after = Asap._N_FACTORIZATIONS[]
        @test after == before                              # no refactorization
        @test sum(dp.c) ≈ dp.δ rtol = 1e-9
        @test sum(dp2.c) ≈ dp2.δ rtol = 1e-9

        # Querying a restrained DOF is a clear error.
        @test_throws ErrorException drift_participation(model, n1.globalID[1])
    end

    # Shared 2-story cantilever column (two stacked elements) for the
    # differential- and batched-virtual-load tests.
    function _two_story_column()
        h = 3.0; E = 200.0e9; I = 8.0e-5; A = 1.0e-2
        sec = Section(A * u"m^2", E * u"Pa", 80.0e9u"Pa",
                      I * u"m^4", I * u"m^4", 1.0e-6u"m^4")
        base = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        mid  = Node([0.0u"m", 0.0u"m", h * 1.0u"m"], :free)
        top  = Node([0.0u"m", 0.0u"m", 2h * 1.0u"m"], :free)
        e1 = Element(base, mid, sec); e2 = Element(mid, top, sec)
        model = FrameModel([base, mid, top], [e1, e2],
            [NodeForce(mid, [5.0e3u"N", 0.0u"N", 0.0u"N"]),
             NodeForce(top, [5.0e3u"N", 0.0u"N", 0.0u"N"])])
        solve!(model)
        return model, mid, top
    end

    @testset "inter-story drift: differential DOF decomposes exactly" begin
        model, mid, top = _two_story_column()
        dtop = top.globalID[1]; dmid = mid.globalID[1]
        # Virtual load +1 at the upper floor, −1 at the lower ⇒ the 2nd-story drift.
        dp = drift_participation(model, dtop, dmid)
        @test dp.δ ≈ model.u[dtop] - model.u[dmid] rtol = 1e-9
        @test sum(dp.c) ≈ dp.δ rtol = 1e-9
        @test abs(dp.nonframe) < 1e-8 * max(1.0, abs(dp.δ))
        # The arbitrary-vector form agrees with the differential overload.
        v = zeros(model.nDOFs); v[dtop] = 1.0; v[dmid] = -1.0
        @test drift_participation(model, v).c ≈ dp.c rtol = 1e-12
    end

    @testset "batched virtual loads match single solves" begin
        model, mid, top = _two_story_column()
        dtop = top.globalID[1]; dmid = mid.globalID[1]
        V = zeros(model.nDOFs, 2)
        V[dtop, 1] = 1.0            # absolute drift at the top
        V[dtop, 2] = 1.0; V[dmid, 2] = -1.0   # 2nd-story inter-story drift
        bdp = drift_participation(model, V)
        @test bdp isa BatchedDriftParticipation
        d1 = drift_participation(model, dtop)
        d2 = drift_participation(model, dtop, dmid)
        @test bdp.c[:, 1] ≈ d1.c rtol = 1e-10
        @test bdp.c[:, 2] ≈ d2.c rtol = 1e-10
        @test bdp.δ[1] ≈ d1.δ rtol = 1e-12
        @test bdp.δ[2] ≈ d2.δ rtol = 1e-12
        @test bdp.element_ids == d1.element_ids
    end

    @testset "rigid diaphragm: participation through the reduced factor" begin
        # Two columns tied to a diaphragm master (the rigid-link reduced solve).
        # The master's X DOF is the story drift; the virtual load must be solved
        # on the SAME reduced system, not the unreduced factor.
        h = 4.0; E = 200.0e9; I = 8.0e-5; A = 1.0e-2; b = 1.5
        sec = Section(A * u"m^2", E * u"Pa", 80.0e9u"Pa",
                      I * u"m^4", I * u"m^4", 1.0e-6u"m^4")
        base1 = Node([0.0u"m",  b * 1.0u"m", 0.0u"m"], :fixed)
        base2 = Node([0.0u"m", -b * 1.0u"m", 0.0u"m"], :fixed)
        top1  = Node([0.0u"m",  b * 1.0u"m", h * 1.0u"m"], :free)
        top2  = Node([0.0u"m", -b * 1.0u"m", h * 1.0u"m"], :free)
        master = Node([0.0u"m", 0.0u"m", h * 1.0u"m"], :free)
        c1 = Element(base1, top1, sec); c2 = Element(base2, top2, sec)
        model = FrameModel([base1, base2, top1, top2, master], [c1, c2],
            [NodeForce(master, [1.0e4u"N", 0.0u"N", 0.0u"N"])])
        links = [RigidLink(master, top1), RigidLink(master, top2)]
        solve!(model, links)

        md = master.globalID[1]
        dp = drift_participation(model, md; links = links)
        @test dp.δ ≈ model.u[md] rtol = 1e-9
        # The links store no strain energy ⇒ the two columns carry all of δ.
        @test sum(dp.c) ≈ dp.δ rtol = 1e-6
        @test abs(dp.nonframe) < 1e-6 * max(1.0, abs(dp.δ))
        @test length(dp.c) == 2
        @test dp.c[1] ≈ dp.c[2] rtol = 1e-6         # symmetric columns participate equally
        @test all(c -> c > 0, dp.c)

        # linear_response on the reduced system reproduces the real displacement.
        @test linear_response(model, model.P .- model.Pf; links = links) ≈ model.u rtol = 1e-7
    end

    @testset "linear_response reproduces the real displacement" begin
        model, _, _ = _two_story_column()
        @test linear_response(model, model.P .- model.Pf) ≈ model.u rtol = 1e-9
        @test_throws ErrorException linear_response(model, zeros(model.nDOFs + 1))
    end

    @testset "linear_response predicts a resize to first order" begin
        # Two identical 2-story columns differing only in element 1's inertia.
        mk(I1) = begin
            h = 3.0; E = 200.0e9; I = 8.0e-5; A = 1.0e-2
            s(Ix) = Section(A * u"m^2", E * u"Pa", 80.0e9u"Pa",
                            Ix * u"m^4", Ix * u"m^4", 1.0e-6u"m^4")
            base = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            mid  = Node([0.0u"m", 0.0u"m", h * 1.0u"m"], :free)
            top  = Node([0.0u"m", 0.0u"m", 2h * 1.0u"m"], :free)
            e1 = Element(base, mid, s(I1)); e2 = Element(mid, top, s(I))
            m = FrameModel([base, mid, top], [e1, e2],
                [NodeForce(mid, [5.0e3u"N", 0.0u"N", 0.0u"N"]),
                 NodeForce(top, [5.0e3u"N", 0.0u"N", 0.0u"N"])])
            solve!(m); m
        end
        I0 = 8.0e-5; ε = 1.0e-3
        m0 = mk(I0); m1 = mk(I0 * (1 + ε))                  # exact perturbed solve
        ΔK = m1.elements[1].K - m0.elements[1].K            # exact element stiffness delta
        g = m0.elements[1].globalID
        f = zeros(m0.nDOFs); f[g] .-= ΔK * m0.u[g]          # perturbation load −ΔK u
        Δu = linear_response(m0, f)                          # first-order prediction
        exact = m1.u .- m0.u
        # The first-order error (O(ε²)) is a small fraction of the step (O(ε)).
        @test norm(exact) > 0
        @test norm((m0.u .+ Δu) .- m1.u) < 1.0e-2 * norm(exact)
    end
end
