#=
Tests for Released-Element Local Stiffness Matrices
===================================================

Validates the member end-release stiffness functions added to
`src/Elements/K.jl`:

    k_fixedfixed, k_fixedfree, k_freefixed, k_freefree, k_joist

These matrices are obtained by static condensation (member end-release
modification) of the fully fixed-fixed 12x12 local stiffness matrix — a
standard closed-form matrix-analysis operation.

References for the known-answer values:
  - W. McGuire, R.H. Gallagher, R.D. Ziemian, "Matrix Structural Analysis,"
    2nd ed., Sec. 7.4 (static condensation).
  - A. Kassimali, "Structural Analysis," 5th ed., Table 8.1 / Table 8.2
    (member stiffness with end moment releases).

Local DOF ordering (node i then node j):
  1: u_x,i  2: u_y,i  3: u_z,i  4: θx,i  5: θy,i  6: θz,i
  7: u_x,j  8: u_y,j  9: u_z,j 10: θx,j 11: θy,j 12: θz,j
Plane assignment (per natural formulation): I3 = Ix ↔ (u_y, θz);  I2 = Iy ↔ (u_z, θy).
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "Released Element Stiffness" begin

    # Independent (test-local) static condensation — deliberately NOT the
    # implementation's helper, so it provides a real cross-check.
    function condense(K::AbstractMatrix, released::Vector{Int})
        n = size(K, 1)
        retained = setdiff(1:n, released)
        Krr = K[retained, retained]
        Krc = K[retained, released]
        Kcc = K[released, released]
        Kred = Krr - Krc * (Kcc \ Krc')
        out = zeros(n, n)
        out[retained, retained] = Kred
        return out
    end

    # Independent fixed-fixed assembly for cross-checks, matching the Asap
    # convention I3 = Ix, I2 = Iy (Bernoulli-Euler: Ay = Az = Inf).
    function full_ff(E, A, L, G, Ix, Iy, J)
        K = zeros(12, 12)
        aN = zeros(6, 12)
        DN = zeros(6, 6)
        Asap.local_stiffness_natural!(K, E, G, A, Iy, Ix, J, Inf, Inf, L, aN, DN)
        return 0.5 * (K + K')
    end

    zero_rows_cols(K, dofs) = all(iszero, K[dofs, :]) && all(iszero, K[:, dofs])

    # (name, E, A, L, G, Ix, Iy, J): a unit case for clean known answers and a
    # realistic steel case with Ix ≠ Iy to catch bending-plane swaps.
    cases = [
        ("unit",      1.0,   1.0,  1.0, 1.0, 1.0,   1.0,   1.0),
        ("steel W",   200e9, 0.01, 3.5, 77e9, 1.0e-4, 5.0e-5, 1.0e-5),
    ]

    for (name, E, A, L, G, Ix, Iy, J) in cases
        EA = E * A / L
        GJ = G * J / L

        @testset "k_fixedfixed [$name]" begin
            K = Asap.k_fixedfixed(E, A, L, G, Ix, Iy, J)
            @test size(K) == (12, 12)
            @test K ≈ K'                                   # symmetry
            # Textbook fixed-fixed coefficients (Kassimali 5e, Table 8.1).
            @test K[1, 1]  ≈ EA
            @test K[1, 7]  ≈ -EA
            @test K[4, 4]  ≈ GJ
            @test K[4, 10] ≈ -GJ
            # Strong plane (I3 = Ix): u_y (2,8) and θz (6,12).
            @test K[2, 2]  ≈ 12 * E * Ix / L^3
            @test K[6, 6]  ≈ 4  * E * Ix / L
            @test K[6, 12] ≈ 2  * E * Ix / L
            @test abs(K[2, 6]) ≈ 6 * E * Ix / L^2
            # Weak plane (I2 = Iy): u_z (3,9) and θy (5,11).
            @test K[3, 3]  ≈ 12 * E * Iy / L^3
            @test K[5, 5]  ≈ 4  * E * Iy / L
            @test abs(K[3, 5]) ≈ 6 * E * Iy / L^2
            # Free-free element: 6 rigid-body modes ⇒ 6 zero eigenvalues.
            ev = sort(eigvals(Symmetric(K)))
            @test count(x -> abs(x) < 1e-6 * (EA + GJ + E*Ix/L^3 + 1), ev) == 6
            @test all(ev .>= -1e-6 * (EA + 1))          # positive semidefinite
            @test K ≈ full_ff(E, A, L, G, Ix, Iy, J)     # matches natural form
        end

        @testset "k_joist [$name]" begin
            K = Asap.k_joist(E, A, L, G, J, Ix, Iy)
            @test size(K) == (12, 12)
            @test K ≈ K'
            # Retained: axial EA/L and torsion GJ/L only (rank 2).
            @test K[1, 1]  ≈ EA
            @test K[1, 7]  ≈ -EA
            @test K[4, 4]  ≈ GJ
            @test K[4, 10] ≈ -GJ
            @test rank(K, atol = 1e-8 * (EA + GJ)) == 2
            # Every bending DOF (translations + rotations, both ends) is zero.
            @test zero_rows_cols(K, [2, 3, 5, 6, 8, 9, 11, 12])
            # Independent check: condensing the 4 bending rotations from the
            # full fixed-fixed matrix collapses bending to zero, leaving
            # axial + torsion — must equal the closed-form k_joist.
            Kc = condense(full_ff(E, A, L, G, Ix, Iy, J), [5, 6, 11, 12])
            @test K ≈ Kc
        end

        @testset "k_freefree [$name]" begin
            K = Asap.k_freefree(E, A, L, Ix, Iy)
            @test size(K) == (12, 12)
            @test K ≈ K'
            # Two-force member: only axial EA/L survives (rank 1).
            @test K[1, 1] ≈ EA
            @test K[1, 7] ≈ -EA
            @test K[7, 7] ≈ EA
            @test rank(K, atol = 1e-8 * EA) == 1
            @test zero_rows_cols(K, [2, 3, 4, 5, 6, 8, 9, 10, 11, 12])
            # Embedded axial block equals Asap's truss member scaling EA/L·[1 -1; -1 1].
            @test [K[1,1] K[1,7]; K[7,1] K[7,7]] ≈ EA * [1.0 -1.0; -1.0 1.0]
        end

        @testset "k_fixedfree [$name]" begin
            K = Asap.k_fixedfree(E, A, L, Ix, Iy)
            @test size(K) == (12, 12)
            @test K ≈ K'
            @test K[1, 1] ≈ EA                             # axial retained
            @test K[1, 7] ≈ -EA
            # All rotations released at end node j ⇒ zero rows/cols there,
            # and torsion released ⇒ both torsion DOFs zero member-wide.
            @test zero_rows_cols(K, [4, 10, 11, 12])
            # Fixed-pinned (propped cantilever), Kassimali 5e Table 8.2:
            # transverse 3EI/L^3, rotational at fixed end 3EI/L.
            @test K[2, 2] ≈ 3 * E * Ix / L^3               # u_y strong plane
            @test K[8, 8] ≈ 3 * E * Ix / L^3
            @test K[2, 8] ≈ -3 * E * Ix / L^3
            @test K[6, 6] ≈ 3 * E * Ix / L                 # θz at fixed end i
            @test K[3, 3] ≈ 3 * E * Iy / L^3               # u_z weak plane
            @test K[5, 5] ≈ 3 * E * Iy / L                 # θy at fixed end i
            # Independent check: condense the full end-j rotation set {10,11,12}
            # (with real G, J) — torsion condenses to zero either way.
            Kc = condense(full_ff(E, A, L, G, Ix, Iy, J), [10, 11, 12])
            @test K ≈ Kc
        end

        @testset "k_freefixed [$name]" begin
            K = Asap.k_freefixed(E, A, L, Ix, Iy)
            @test size(K) == (12, 12)
            @test K ≈ K'
            @test K[1, 1] ≈ EA
            @test K[1, 7] ≈ -EA
            # All rotations released at start node i; torsion zero member-wide.
            @test zero_rows_cols(K, [4, 5, 6, 10])
            # Pinned-fixed (mirror of fixed-free): rotational stiffness now at
            # the fixed end j.
            @test K[2, 2]   ≈ 3 * E * Ix / L^3
            @test K[8, 8]   ≈ 3 * E * Ix / L^3
            @test K[12, 12] ≈ 3 * E * Ix / L               # θz at fixed end j
            @test K[3, 3]   ≈ 3 * E * Iy / L^3
            @test K[11, 11] ≈ 3 * E * Iy / L               # θy at fixed end j
            Kc = condense(full_ff(E, A, L, G, Ix, Iy, J), [4, 5, 6])
            @test K ≈ Kc
        end
    end

    # =========================================================================
    # Dispatch / units integration: local_K must resolve for every release type
    # (this is the seam that previously threw UndefVarError: k_joist).
    # =========================================================================
    @testset "local_K dispatch for all release types" begin
        E = 200e9u"Pa"; G = 77e9u"Pa"
        A = 0.01u"m^2"; Ix = 1e-4u"m^4"; Iy = 5e-5u"m^4"; J = 1e-5u"m^4"
        section = Section(A, E, G, Ix, Iy, J)
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([3.5u"m", 0.0u"m", 0.0u"m"], :free)

        for rel in (:fixedfixed, :fixedfree, :freefixed, :freefree, :joist)
            elem = Element(n1, n2, section; release = rel)
            elem.length = 3.5u"m"   # populated by process! in a real model
            K = Asap.local_K(elem)
            @test size(K) == (12, 12)
            @test all(isfinite, K)
            @test K ≈ K'
        end

        # Joist explicitly retains axial + torsion in SI units (EA/L, GJ/L).
        joist = Element(n1, n2, section; release = :joist)
        joist.length = 3.5u"m"
        Kj = Asap.local_K(joist)
        L = 3.5
        @test Kj[1, 1]  ≈ 200e9 * 0.01 / L
        @test Kj[4, 4]  ≈ 77e9 * 1e-5 / L
        @test all(iszero, Kj[[2, 3, 5, 6, 8, 9, 11, 12], :])
    end
end
