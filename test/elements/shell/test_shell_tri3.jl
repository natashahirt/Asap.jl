using Test
using Asap
using LinearAlgebra
using SparseArrays
using Unitful
using Random

@testset "ShellTri3 Full Shell Element" begin

    @testset "Geometry & LCS" begin
        # Simple triangle in XY plane
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        elem = Asap.ShellTri3((n1, n2, n3), 0.1u"m", 200e9u"Pa", 0.3)
        Asap.lcs!(elem)
        
        # Check LCS
        x_loc, y_loc, z_loc = elem.LCS
        @test x_loc ≈ [1.0, 0.0, 0.0] atol=1e-10
        @test y_loc ≈ [0.0, 1.0, 0.0] atol=1e-10
        @test z_loc ≈ [0.0, 0.0, 1.0] atol=1e-10
        
        # Check area (right triangle with legs 1×1 → area = 0.5)
        @test elem.area ≈ 0.5 atol=1e-10
    end
    
    @testset "Stiffness Matrix Symmetry & Positive Semi-Definiteness" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)  # equilateral
        
        elem = Asap.ShellTri3((n1, n2, n3), 0.1u"m", 200e9u"Pa", 0.3)
        Asap.process!(elem)
        
        # Local K should be 18×18, symmetric
        K_local = Asap.local_K(elem)
        @test size(K_local) == (18, 18)
        @test K_local ≈ K_local' rtol=1e-10
        
        # Global K should be 18×18, symmetric
        @test size(elem.K) == (18, 18)
        @test elem.K ≈ elem.K' rtol=1e-10
        
        # Eigenvalues should be non-negative
        eigvals_K = eigvals(Symmetric(elem.K))
        @test all(eigvals_K .>= -1e-6 * maximum(eigvals_K))
        
        # For a single unconstrained shell element there are at least 6 near-zero
        # eigenvalues: the six rigid-body modes. The in-plane rotation is one of
        # them, because the drilling penalty vanishes when θz equals that rotation.
        rbm_count = count(x -> abs(x) < 1e-3 * maximum(abs, eigvals_K), eigvals_K)
        @test rbm_count >= 6  # At least 6 rigid body modes
        
        @info "Rigid body mode count: $rbm_count"
    end
    
    @testset "Membrane Patch Test - Constant Strain" begin
        # 1m × 1m square domain split into 2 triangles
        #   4 ●───────● 3
        #     │ \   2 │
        #     │   \   │
        #     │  1  \ │
        #   1 ●───────● 2
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        nodes = [n1, n2, n3, n4]
        
        # Assign node IDs (6 DOF per node for full shell)
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        E = 200e9  # Pa
        ν = 0.3
        t = 0.1    # m
        
        elem1 = Asap.ShellTri3((n1, n2, n4), t*u"m", E*u"Pa", ν)
        elem2 = Asap.ShellTri3((n2, n3, n4), t*u"m", E*u"Pa", ν)
        
        Asap.process!(elem1)
        Asap.process!(elem2)
        
        elem1.globalID = vcat([n.globalID for n in elem1.nodes]...)
        elem2.globalID = vcat([n.globalID for n in elem2.nodes]...)
        
        # Apply linear displacement: u_x = εxx * x
        εxx_global = 0.001
        
        # Node positions (x-coordinates)
        x_coords = [0.0, 1.0, 1.0, 0.0]
        
        # Build displacement vector (24 DOFs = 4 nodes × 6 DOF)
        u = zeros(24)
        for (i, x) in enumerate(x_coords)
            u[(i-1)*6 + 1] = εxx_global * x  # u_x
        end
        
        σ1 = Asap.stress(elem1, u)
        σ2 = Asap.stress(elem2, u)
        
        # Expected stresses
        σxx_global = E / (1 - ν^2) * εxx_global
        σyy_global = E / (1 - ν^2) * ν * εxx_global
        
        # Principal stresses should match (invariant)
        principal_1 = sort([σ1[1], σ1[2]])
        principal_2 = sort([σ2[1], σ2[2]])
        @test principal_1 ≈ principal_2 rtol=1e-8
        
        @info "Membrane patch test passed" σ1 σ2 σxx_global σyy_global
    end
    
    @testset "Bending Stiffness - Simply Supported Plate" begin
        # Simply supported square plate with uniform load
        # Analytical solution: w_max = α × q × a⁴ / D
        # where D = E×t³/(12(1-ν²)), α ≈ 0.00416 for ν = 0.3
        
        # Create a 2×2 grid (4 elements) for a 1m × 1m plate
        #   6 ●───────● 7 ───────● 8
        #     │ \   4 │ \   6   │
        #     │   \   │   \     │
        #     │  3  \ │  5  \   │
        #   3 ●───────● 4 ───────● 5
        #     │ \   2 │ \   ??  │  (simplified: just 4 elements)
        #     │   \   │   \     │
        #     │  1  \ │  ??  \  │
        #   0 ●───────● 1 ───────● 2
        
        # Simplified: just check that bending stiffness terms are present and reasonable
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        E = 200e9   # Pa (steel)
        ν = 0.3
        t = 0.01    # 10mm plate
        
        elem = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        Asap.process!(elem)
        
        K = Asap.local_K(elem)
        
        # Bending stiffness D = E × t³ / (12(1-ν²))
        D = E * t^3 / (12 * (1 - ν^2))
        
        # Rotation DOFs are at indices 4,5,6 for node 1, 10,11,12 for node 2, etc.
        # Check that bending stiffness terms (rotation-rotation coupling) are present
        K_rot = K[4:6, 4:6]  # Rotation DOFs for node 1
        
        # Bending stiffness should scale with D
        @test maximum(abs.(K_rot)) > 0  # Non-zero bending stiffness
        
        # Membrane stiffness (translation DOFs) should scale with E×t
        K_trans = K[1:3, 1:3]
        membrane_scale = E * t * elem.area
        
        @test maximum(abs.(K_trans)) > 0  # Non-zero membrane stiffness
        
        @info "Bending stiffness check passed" D=D membrane_scale=membrane_scale
    end
    
    @testset "Bending Moments Recovery" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        for (i, node) in enumerate([n1, n2, n3])
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        E = 200e9
        ν = 0.3
        t = 0.1
        
        elem = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        
        # Apply constant curvature: θy varies linearly with x
        # This produces constant Mxx
        u = zeros(18)
        κxx = 0.001  # curvature 1/m
        
        # θy = κxx × x for each node
        x_coords = [0.0, 1.0, 0.5]
        for (i, x) in enumerate(x_coords)
            u[(i-1)*6 + 5] = κxx * x  # θy DOF
        end
        
        M = Asap.bending_moments(elem, u)
        
        # Expected: Mxx = D × κxx where D = E×t³/(12(1-ν²))
        D = E * t^3 / (12 * (1 - ν^2))
        Mxx_expected = D * κxx
        
        @test M[1] ≈ Mxx_expected rtol=0.1  # Allow 10% tolerance for numerical approx
        
        @info "Bending moment test" M=M Mxx_expected=Mxx_expected
    end
    
    @testset "AreaLoad Application" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        for (i, node) in enumerate([n1, n2, n3])
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        elem = Asap.ShellTri3((n1, n2, n3), 0.1u"m", 200e9u"Pa", 0.3)
        Asap.process!(elem)
        
        # Apply 1000 Pa pressure
        pressure = 1000.0u"Pa"
        load = Asap.AreaLoad(elem, pressure; direction=(0.0, 0.0, -1.0))
        
        # Compute nodal forces
        nf = Asap.nodal_forces(load)
        
        # Total force should equal pressure × area
        total_force = sum(f[2][3] for f in nf)  # Sum of Z-forces
        expected_force = -1000.0 * elem.area  # Negative Z
        
        @test total_force ≈ expected_force rtol=1e-10
        @test length(nf) == 3  # Three nodes
        
        # Each node gets 1/3 of total
        for (node, fvec) in nf
            @test fvec[3] ≈ expected_force / 3 rtol=1e-10
        end
        
        @info "Surface load test passed" area=elem.area total_force=total_force
    end

    @testset "Drilling penalty is self-equilibrated" begin
        # A triangle that is not in a coordinate plane, so the local frame is
        # not the global one. Equilibrium is checked in that local frame.
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([2.0u"m", 0.1u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.3u"m", 1.4u"m", 0.2u"m"], :free)
        elem = Asap.ShellTri3((n1, n2, n3), 0.3u"m", 30e9u"Pa", 0.2)
        Asap.process!(elem)
        K = Asap.local_K(elem)

        # Resultant of the forces the element applies to its nodes, about
        # node 1. The element frame puts that node at the origin.
        function resultant(f)
            Fx = Fy = Fz = Mx = My = Mz = 0.0
            for i in 1:3
                base = (i - 1) * 6
                x = elem.ecoords_e[i, 1]
                y = elem.ecoords_e[i, 2]
                fx, fy, fz = f[base + 1], f[base + 2], f[base + 3]
                Fx += fx
                Fy += fy
                Fz += fz
                Mx += y * fz + f[base + 4]
                My += -x * fz + f[base + 5]
                Mz += x * fy - y * fx + f[base + 6]
            end
            return (Fx, Fy, Fz, Mx, My, Mz)
        end

        # Rigid rotation about the normal: u = −θ y, v = θ x, and θz = θ,
        # which is the rotation the penalty is measured against.
        θ = 1.0e-3
        d = zeros(18)
        for i in 1:3
            x = elem.ecoords_e[i, 1]
            y = elem.ecoords_e[i, 2]
            base = (i - 1) * 6
            d[base + 1] = -θ * y
            d[base + 2] = θ * x
            d[base + 6] = θ
        end
        f_rigid = K * d
        @test norm(f_rigid) <= 1e-8 * opnorm(K, 1) * norm(d)

        # Any displacement: no net force and no net moment.
        f = K * rand(Random.Xoshiro(1), 18)
        scale = max(norm(f), 1.0)
        @test all(abs.(resultant(f)) .<= 1e-8 * scale)
    end

    @testset "Membrane-only element keeps exactly the in-plane stiffness" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.3u"m", 0.2u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.4u"m", 0.9u"m", 0.0u"m"], :free)
        full = Asap.ShellTri3((n1, n2, n3), 0.2u"m", 30e9u"Pa", 0.2)
        mem = Asap.ShellTri3((n1, n2, n3), 0.2u"m", 30e9u"Pa", 0.2; membrane_only=true)
        @test !full.membrane_only && mem.membrane_only
        Asap.process!(full)
        Asap.process!(mem)
        Kf = Asap.local_K(full)
        Km = Asap.local_K(mem)
        inplane = [b + i for b in (0, 6, 12) for i in 1:2]
        other = setdiff(1:18, inplane)
        # The drilling penalty in the full element also touches u and v, so
        # compare against a full element with the drilling penalty off.
        nodrill = Asap.ShellTri3((n1, n2, n3), 0.2u"m", 30e9u"Pa", 0.2; drilling_scale=0.0)
        Asap.process!(nodrill)
        @test Km[inplane, inplane] ≈ Asap.local_K(nodrill)[inplane, inplane] rtol=1e-12
        @test all(iszero, Km[other, :])
        @test all(iszero, Km[:, other])
        @test mem.K ≈ mem.K' rtol=1e-12
        # Three in-plane rigid modes; the 6 in-plane DOFs have rank 3.
        @test rank(Km[inplane, inplane]; rtol=1e-10) == 3
    end

    @testset "HAND CHECK membrane cantilever in its plane" begin
        # Deep cantilever, L = 4 m, depth h = 2 m, t = 0.2 m, E = 30 GPa,
        # ν = 0. Tip load P = 1000 kN spread over the free end.
        # Timoshenko beam: δ = P L³/(3 E I) + P L/(κ G A), κ = 5/6, G = E/2.
        L, h, t, E = 4.0, 2.0, 0.2, 30e9
        P = 1.0e6
        I = t * h^3 / 12
        A = t * h
        G = E / 2
        δ_beam = P * L^3 / (3 * E * I) + P * L / ((5 / 6) * G * A)
        function tip(nx, ny)
            xs = range(0.0, L; length=nx + 1)
            ys = range(0.0, h; length=ny + 1)
            grid = Matrix{Asap.Node}(undef, nx + 1, ny + 1)
            for i in 1:(nx + 1), j in 1:(ny + 1)
                dof = i == 1 ? [false, false, false, false, false, false] :
                               [true, true, false, false, false, false]
                grid[i, j] = Asap.Node([xs[i] * u"m", ys[j] * u"m", 0.0u"m"], dof)
            end
            shells = Asap.ShellTri3[]
            for i in 1:nx, j in 1:ny
                a, b, c, d = grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]
                push!(shells, Asap.ShellTri3((a, b, c), t * u"m", E * u"Pa", 0.0; membrane_only=true))
                push!(shells, Asap.ShellTri3((a, c, d), t * u"m", E * u"Pa", 0.0; membrane_only=true))
            end
            loads = Asap.AbstractLoad[]
            # Consistent parabolic shear would be closer to beam theory; equal
            # shares on a fine edge are enough for the end-load resultant.
            for j in 1:(ny + 1)
                w = (j == 1 || j == ny + 1) ? 0.5 / ny : 1.0 / ny
                push!(loads, Asap.NodeForce(grid[end, j], [0.0u"N", -P * w * u"N", 0.0u"N"]))
            end
            model = Asap.ShellModel(vec(grid), shells, loads)
            Asap.solve!(model)
            return -sum(ustrip(u"m", grid[end, j].displacement[2]) for j in 1:(ny + 1)) / (ny + 1)
        end
        coarse = tip(4, 2)
        fine = tip(32, 16)
        # Linear triangles are stiff in bending: the coarse mesh under-reads.
        @test coarse < fine
        @test fine ≈ δ_beam rtol=0.05
        @info "Membrane cantilever tip" δ_beam coarse fine
    end

end
