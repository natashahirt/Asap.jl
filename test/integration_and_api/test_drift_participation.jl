using Asap
using Test
using Unitful
using LinearAlgebra: dot

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
end
