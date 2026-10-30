#=
Tests for the RigidLink kinematic-offset element
================================================

RigidLink ties a slave node to a master node through a rigid body:
    u_B = u_A + θ_A × r_AB ,  θ_B = θ_A
resolved by master/slave DOF reduction in `solve!(model, links)`
(Bathe, "Finite Element Procedures" §4.4; Cook et al. §6.3).

Every case below is checked against a CLOSED-FORM analytical result. All beams
use the default (Bernoulli-Euler, shear-rigid) `Section`, so the classical
cantilever relations are exact and tolerances are tight.

Common cantilever relations used throughout (fixed base, height h, bending EI):
- Tip lateral flexibility under tip shear F and tip moment M:
      v = F h³/(3EI) + M h²/(2EI)
      φ = F h²/(2EI) + M h/(EI)
- Tip 2×2 stiffness (invert the above):
      [F; M] = (EI) * [[12/h³, -6/h²], [-6/h², 4/h]] * [v; φ]
=#

using Test
using Asap
using Unitful
using LinearAlgebra
using SparseArrays: spzeros

@testset "RigidLink" begin

    # ---- shared material / section (symmetric so bending axis is irrelevant) ----
    E   = 200e9u"Pa"
    Gm  = E / 2.6
    A   = 0.01u"m^2"
    Iv  = 1e-4u"m^4"          # Ix = Iy = I  → same EI about either bending axis
    J   = 1e-5u"m^4"
    sec = Section(A, E, Gm, Iv, Iv, J)

    E_val  = 200e9            # Pa
    I_val  = 1e-4             # m^4
    EI     = E_val * I_val    # = 2e7 N·m²
    h      = 4.0              # column height [m]

    # =========================================================================
    # Sanity: the kinematic transformation matrix itself
    # =========================================================================
    @testset "Kinematic transformation matrix" begin
        # Master at origin, slave offset r_AB = (2, 3, 4).
        nA = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        nB = Node([2.0u"m", 3.0u"m", 4.0u"m"], :free)
        T = rigid_link_transformation(nA, nB)

        # Rotation rows are identity: θ_B = θ_A
        @test T[4:6, 4:6] ≈ I(3)
        # Translation ← translation is identity
        @test T[1:3, 1:3] ≈ I(3)
        # Translation ← master rotation is θ_A × r_AB  (r = (rx,ry,rz)=(2,3,4)):
        #   u_Bx = ... + θy·rz − θz·ry ;  u_By = ... + θz·rx − θx·rz ;  u_Bz = ... + θx·ry − θy·rx
        @test T[1, 5] ≈ 4.0    # rz
        @test T[1, 6] ≈ -3.0   # -ry
        @test T[2, 4] ≈ -4.0   # -rz
        @test T[2, 6] ≈ 2.0    # rx
        @test T[3, 4] ≈ 3.0    # ry
        @test T[3, 5] ≈ -2.0   # -rx

        # Rigid-body check: a pure master rotation θ about z moves the slave by θ×r.
        θ = [0.0, 0.0, 0.1]                      # rad about global z
        uθ = T * [0.0, 0.0, 0.0, θ...]
        r  = [2.0, 3.0, 4.0]
        @test uθ[1:3] ≈ cross(θ, r)              # translation = θ × r
        @test uθ[4:6] ≈ θ                        # rotation unchanged
    end

    # =========================================================================
    # (a) Rigid diaphragm under a lateral point load
    # -------------------------------------------------------------------------
    # N identical cantilever columns (height h, bending EI) stand along the
    # global y-axis; their tops are slaved (fully rigid) to a master node at the
    # story centroid. A lateral load H (global x) is applied at the master.
    #
    # Derivation: columns share the master's lateral translation v and rotation
    # φ (about y). With the master free to rotate and no external moment, the
    # shared-rotation equilibrium (Σ tip-moment = 0) reduces each column to a
    # FREE-TIP cantilever (3EI/h³). N in parallel ⇒ K = 3N·EI/h³, hence
    #       v = H h³ / (3 N EI) ,   φ = 1.5 v / h .
    # Because the columns lie on the y-axis (rx = 0) the master's φ produces NO
    # vertical (axial) motion of the tops (u_z = −φ·rx = 0), so there is no
    # axial coupling and the free-tip result is exact.
    # =========================================================================
    @testset "(a) Rigid diaphragm — lateral point load" begin
        H   = 10_000.0                 # N
        b   = 1.5                      # column half-spacing along y [m]
        N   = 2

        base1 = Node([0.0u"m",  b*u"m", 0.0u"m"], :fixed)
        base2 = Node([0.0u"m", -b*u"m", 0.0u"m"], :fixed)
        top1  = Node([0.0u"m",  b*u"m", h*u"m"], :free)
        top2  = Node([0.0u"m", -b*u"m", h*u"m"], :free)
        master = Node([0.0u"m", 0.0u"m", h*u"m"], :free)

        nodes = [base1, base2, top1, top2, master]
        elements = [Element(base1, top1, sec), Element(base2, top2, sec)]
        loads = [NodeForce(master, [H*u"N", 0.0u"N", 0.0u"N"])]

        model = FrameModel(nodes, elements, loads)
        links = [RigidLink(master, top1), RigidLink(master, top2)]
        solve!(model, links)

        v_expected = H * h^3 / (3 * N * EI)      # = 10000*64/(3*2*2e7) = 5.3333e-3 m
        φ_expected = 1.5 * v_expected / h

        ux_master = model.u[master.globalID[1]]
        θy_master = model.u[master.globalID[5]]

        @test isapprox(ux_master, v_expected; rtol = 1e-6)
        @test isapprox(abs(θy_master), φ_expected; rtol = 1e-6)

        # Rigid diaphragm: every slaved top translates identically with the master
        # (rx = 0 ⇒ no twist contribution) and does NOT move vertically.
        for top in (top1, top2)
            @test isapprox(model.u[top.globalID[1]], ux_master; rtol = 1e-9)
            @test abs(model.u[top.globalID[3]]) < 1e-10       # u_z ≈ 0 (no axial coupling)
        end
    end

    # =========================================================================
    # (b) Eccentric brace working point
    # -------------------------------------------------------------------------
    # A single vertical cantilever column; its top C is the master. A brace
    # working point W is offset e in +x from C (same height). A vertical force P
    # (global −z, i.e. axial to the column) is applied at W.
    #
    # The rigid collar transfers the eccentric force to the column node as an
    # axial force P plus a moment r_CW × F about y:  M = e·P.  A pure tip moment
    # on a cantilever is constant along its length, so the base moment about y is
    #       M_base = e · P            (base shear = 0, base axial = P).
    # =========================================================================
    @testset "(b) Eccentric brace — F×e moment transfer" begin
        P = 20_000.0                   # N (magnitude of vertical force)
        e = 0.5                        # eccentricity in x [m]

        base = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        C    = Node([0.0u"m", 0.0u"m", h*u"m"], :free)        # master (column top)
        W    = Node([e*u"m", 0.0u"m", h*u"m"], :free)          # slave (working point)

        nodes = [base, C, W]
        elements = [Element(base, C, sec)]
        loads = [NodeForce(W, [0.0u"N", 0.0u"N", -P*u"N"])]

        model = FrameModel(nodes, elements, loads)
        links = [RigidLink(C, W)]
        solve!(model, links)

        M_base_expected = e * P                    # = 10 000 N·m

        # Base reactions (fixed node): moment about y, vertical force, horizontal shear.
        My_base = model.reactions[base.globalID[5]]
        Fz_base = model.reactions[base.globalID[3]]
        Fx_base = model.reactions[base.globalID[1]]

        @test isapprox(abs(My_base), M_base_expected; rtol = 1e-6)
        @test isapprox(Fz_base, P; rtol = 1e-6)                # support pushes up by P
        @test abs(Fx_base) < 1e-6 * P                          # no horizontal reaction
    end

    # =========================================================================
    # (c) Vertical slave node with an applied moment
    # -------------------------------------------------------------------------
    # Single cantilever column, top C is the master. Slave S is a vertical
    # offset a above C. A moment M (about y) is applied at S.
    #
    # With rotations tied, the moment transmits fully to the master node, so the
    # column sees a pure tip moment M:
    #       θ_master = M h / (EI) ,   u_master,x = M h² / (2 EI) .
    # Kinematics then fix the slave:  θ_S = θ_master and
    #       u_S,x = u_master,x + θ_master · a   (since r_CS = (0,0,a), rz = a).
    # =========================================================================
    @testset "(c) Vertical slave — applied moment transmits to master" begin
        M = 15_000.0                   # N·m
        a = 1.0                        # vertical offset of slave above master [m]

        base = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        C    = Node([0.0u"m", 0.0u"m", h*u"m"], :free)          # master
        S    = Node([0.0u"m", 0.0u"m", (h+a)*u"m"], :free)      # slave (above)

        nodes = [base, C, S]
        elements = [Element(base, C, sec)]
        loads = [NodeMoment(S, [0.0u"N*m", M*u"N*m", 0.0u"N*m"])]

        model = FrameModel(nodes, elements, loads)
        links = [RigidLink(C, S)]
        solve!(model, links)

        θ_expected  = M * h / EI                   # 0.003 rad
        ux_expected = M * h^2 / (2 * EI)           # 0.006 m
        uxS_expected = ux_expected + θ_expected * a   # 0.009 m

        θy_master = model.u[C.globalID[5]]
        ux_master = model.u[C.globalID[1]]
        θy_slave  = model.u[S.globalID[5]]
        ux_slave  = model.u[S.globalID[1]]

        @test isapprox(abs(θy_master), θ_expected;  rtol = 1e-6)
        @test isapprox(abs(ux_master), ux_expected; rtol = 1e-6)
        @test isapprox(θy_slave, θy_master; rtol = 1e-9)                 # θ_B = θ_A
        @test isapprox(abs(ux_slave), uxS_expected; rtol = 1e-6)        # u_B = u_A + θ_A×r

        # The applied moment reaches the base (pure tip moment ⇒ constant along span).
        @test isapprox(abs(model.reactions[base.globalID[5]]), M; rtol = 1e-6)
    end

    # =========================================================================
    # (d) Per-DOF release — rotational decoupling
    # -------------------------------------------------------------------------
    # Two coincident cantilever columns (r_AB = 0): column A top is the master
    # C, column B top is the slave S. A moment M (about y) is applied at S.
    #
    # FULLY RIGID (θ_y tied): C and S share the rotation; the two columns act in
    # parallel under the top moment M:
    #       θ_C = θ_S = M h / (2 EI) ,  master base moment = M / 2 .
    #
    # θ_y RELEASED (release = (F,F,F,F,T,F)): the slave rotation is independent.
    # Translations stay tied (r = 0 ⇒ shared u_x). Solving the 3-DOF reduced
    # system (u, θ_C, θ_S) with M_A = 0 (no master moment), M_B = M (slave), and
    # F_A + F_B = 0 gives:
    #       u = M h²/(4EI) ,  θ_C = 1.5 u/h ,  θ_S = 2.5 u/h ,
    #       master base moment = 3 M / 4 .
    # The key qualitative result: the released axis DECOUPLES the rotations, so
    # θ_C ≠ θ_S (whereas they are equal in the rigid case).
    # =========================================================================
    @testset "(d) Per-DOF release — rotational decoupling" begin
        M = 15_000.0                   # N·m

        function two_column_model(release)
            baseA = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            baseB = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            C     = Node([0.0u"m", 0.0u"m", h*u"m"], :free)     # master (col A top)
            S     = Node([0.0u"m", 0.0u"m", h*u"m"], :free)     # slave  (col B top)
            nodes = [baseA, baseB, C, S]
            elements = [Element(baseA, C, sec), Element(baseB, S, sec)]
            loads = [NodeMoment(S, [0.0u"N*m", M*u"N*m", 0.0u"N*m"])]
            model = FrameModel(nodes, elements, loads)
            links = [RigidLink(C, S; release = release)]
            return model, links, baseA, C, S
        end

        # --- fully rigid ---
        m_r, l_r, baseA_r, C_r, S_r = two_column_model((false, false, false, false, false, false))
        solve!(m_r, l_r)
        θC_rigid = m_r.u[C_r.globalID[5]]
        θS_rigid = m_r.u[S_r.globalID[5]]
        My_base_rigid = m_r.reactions[baseA_r.globalID[5]]

        @test isapprox(θC_rigid, θS_rigid; rtol = 1e-9)                 # tied ⇒ equal
        @test isapprox(abs(θC_rigid), M * h / (2 * EI); rtol = 1e-6)    # 0.0015 rad
        @test isapprox(abs(My_base_rigid), M / 2; rtol = 1e-6)         # 7 500 N·m

        # --- θ_y released ---
        m_x, l_x, baseA_x, C_x, S_x = two_column_model((false, false, false, false, true, false))
        solve!(m_x, l_x)
        θC_rel = m_x.u[C_x.globalID[5]]
        θS_rel = m_x.u[S_x.globalID[5]]
        My_base_rel = m_x.reactions[baseA_x.globalID[5]]

        u_expected  = M * h^2 / (4 * EI)           # 0.003 m
        θC_expected = 1.5 * u_expected / h         # 0.001125 rad
        θS_expected = 2.5 * u_expected / h         # 0.001875 rad

        # Decoupled: master and slave rotations differ.
        @test !isapprox(θC_rel, θS_rel; rtol = 1e-3)
        @test isapprox(abs(θC_rel), θC_expected; rtol = 1e-6)
        @test isapprox(abs(θS_rel), θS_expected; rtol = 1e-6)
        @test isapprox(abs(My_base_rel), 3 * M / 4; rtol = 1e-6)       # 11 250 N·m

        # The released axis changes the master response relative to the tied case.
        @test !isapprox(abs(θC_rel), abs(θC_rigid); rtol = 1e-3)
    end

    # =========================================================================
    # (e) Story geometric stiffness on the master (multi-case solve)
    # -------------------------------------------------------------------------
    # Four cantilever columns at (±b, ±b) under a rigid floor with the master at
    # the centre. A story weight W rests on the columns, W/4 each.
    #
    # Wilson & Habibullah (1987), Eq. 7: the story couple is (W/h)[1 −1; −1 1]
    # on the story drift. With a fixed base only the top term survives, so
    # K_G[ux, ux] = K_G[uy, uy] = −W/h on the master. Eq. 11 with the exact
    # Σ P_m r_m² form gives K_G[θz, θz] = −Σ (W/4)(2b²)/h = −W·2b²/h.
    #
    # A rank-one change d on one DOF changes that DOF's flexibility f to
    # f/(1 + d f) (Sherman–Morrison), so with u₀ = f·F the second-order answer is
    #       u = u₀ / (1 − θ) ,  θ = |d| u₀ / F ,
    # exactly, whatever the elastic coupling to the other master DOFs.
    # For translation θ = W Δ₀ / (H h): the ASCE 7-22 §12.8.7 stability ratio.
    # =========================================================================
    @testset "(e) Geometric stiffness on the master — 1/(1−θ)" begin
        b = 1.5
        W = 2.0e6                       # N, story weight
        H = 10_000.0                    # N
        Tz = 20_000.0                   # N·m

        corners = [(b, b), (-b, b), (-b, -b), (b, -b)]
        bases = [Node([x*u"m", y*u"m", 0.0u"m"], :fixed) for (x, y) in corners]
        tops = [Node([x*u"m", y*u"m", h*u"m"], :free) for (x, y) in corners]
        master = Node([0.0u"m", 0.0u"m", h*u"m"], :free)
        nodes = vcat(bases, tops, [master])
        elements = [Element(bases[i], tops[i], sec) for i in 1:4]
        model = Model(nodes, elements, AbstractLoad[])
        links = [RigidLink(master, t) for t in tops]
        process!(model)

        cases = [AbstractLoad[NodeForce(master, [H*u"N", 0.0u"N", 0.0u"N"])],
                 AbstractLoad[NodeMoment(master, [0.0u"N*m", 0.0u"N*m", Tz*u"N*m"])]]
        u1 = solve(model, cases, links)

        gx = master.globalID[1]
        gy = master.globalID[2]
        gθ = master.globalID[6]
        dx = -W / h
        dθ = -W * 2b^2 / h
        Kg = spzeros(model.nDOFs, model.nDOFs)
        Kg[gx, gx] = dx
        Kg[gy, gy] = dx
        Kg[gθ, gθ] = dθ
        u2 = solve(model, cases, links; geometric_stiffness = Kg)

        Δ0 = u1[1][gx]
        θ_sway = abs(dx) * Δ0 / H
        @test 0.0 < θ_sway < 0.5
        @test isapprox(u2[1][gx], Δ0 / (1 - θ_sway); rtol = 1e-9)

        φ0 = u1[2][gθ]
        θ_twist = abs(dθ) * φ0 / Tz
        @test 0.0 < θ_twist < 0.5
        @test isapprox(u2[2][gθ], φ0 / (1 - θ_twist); rtol = 1e-9)

        # Slaves still follow the rigid floor in the second-order answer.
        for t in tops
            ry = Asap.to_meters(t.position[2])
            @test isapprox(u2[2][t.globalID[1]], -ry * u2[2][gθ]; rtol = 1e-9)
        end

        # The base shear carries H plus the P-delta shear W Δ / h.
        u = u2[1]
        Vbase = 0.0
        for el in elements
            ue = u[el.globalID]
            f = el.K * ue
            Vbase += f[1]
        end
        @test isapprox(abs(Vbase), H + W * u[gx] / h; rtol = 1e-9)

        # Terms on a tied slave DOF are refused.
        Kbad = spzeros(model.nDOFs, model.nDOFs)
        Kbad[tops[1].globalID[1], tops[1].globalID[1]] = -1.0
        @test_throws ArgumentError solve(model, cases, links; geometric_stiffness = Kbad)
        @test_throws DimensionMismatch solve(model, cases, links;
                                             geometric_stiffness = spzeros(3, 3))
    end

end
