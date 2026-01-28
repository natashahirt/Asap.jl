using Test
using Asap
using LinearAlgebra
using SparseArrays
using Unitful

@testset "ShellTri3 Element" begin

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
        
        # Local K should be 6×6, symmetric (use relative tolerance for large values)
        K_local = Asap.local_K(elem)
        @test size(K_local) == (6, 6)
        @test K_local ≈ K_local' rtol=1e-10
        
        # Global K should be 9×9, symmetric
        @test size(elem.K) == (9, 9)
        @test elem.K ≈ elem.K' rtol=1e-10
        
        # Eigenvalues should be non-negative
        # A membrane in 3D has 6 rigid body modes:
        #   - 3 translations (X, Y, Z global)
        #   - 3 rotations (but membrane has no bending stiffness, so out-of-plane rotation is free)
        # Local 6×6 has rank 3 (6 DOF - 3 RBM in plane)
        # Global 9×9 = R'(6×9) × K_local(6×6) × R(6×9) has rank ≤ 3
        eigvals_K = eigvals(Symmetric(elem.K))
        @test all(eigvals_K .>= -1e-6 * maximum(eigvals_K))  # scaled tolerance
        @test count(x -> abs(x) < 1e-3 * maximum(abs, eigvals_K), eigvals_K) == 6  # 6 near-zero modes
    end
    
    @testset "Patch Test - Constant Strain (2D in XY plane)" begin
        # 1m × 1m square domain split into 2 triangles
        # Test with displacement control to verify constant strain recovery
        #
        #   4 ●───────● 3
        #     │ \   2 │
        #     │   \   │
        #     │  1  \ │
        #   1 ●───────● 2
        #
        # For patch test: apply linear displacement field u = εxx * x
        # This should produce constant strain εxx throughout
        
        # Nodes (in XY plane, z=0)
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        nodes = [n1, n2, n3, n4]
        
        # Assign node IDs
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*3+1 : i*3)  # 3 DOF per node for translations
        end
        
        # Material
        E = 200e9  # Pa
        ν = 0.3
        t = 0.1    # m
        
        # Two triangles covering the square
        elem1 = Asap.ShellTri3((n1, n2, n4), t*u"m", E*u"Pa", ν)
        elem2 = Asap.ShellTri3((n2, n3, n4), t*u"m", E*u"Pa", ν)
        
        # Process elements
        Asap.process!(elem1)
        Asap.process!(elem2)
        
        # Assign element global IDs
        elem1.globalID = vcat(n1.globalID, n2.globalID, n4.globalID)
        elem2.globalID = vcat(n2.globalID, n3.globalID, n4.globalID)
        
        # Apply linear displacement: u_x = εxx * x, u_y = 0, u_z = 0
        εxx_global = 0.001  # 0.1% strain in global X direction
        
        # Node positions (x-coordinates)
        x1, x2, x3, x4 = 0.0, 1.0, 1.0, 0.0
        
        # Prescribed displacements
        u = zeros(12)
        u[1] = εxx_global * x1  # node 1 x
        u[4] = εxx_global * x2  # node 2 x  
        u[7] = εxx_global * x3  # node 3 x
        u[10] = εxx_global * x4 # node 4 x
        # All y, z displacements are zero
        
        # Compute stresses (in LOCAL element coordinates)
        σ1 = Asap.stress(elem1, u)
        σ2 = Asap.stress(elem2, u)
        
        # Expected: plane stress in GLOBAL coords: εxx = 0.001, εyy = 0, γxy = 0
        # σxx_global = E/(1-ν²) * εxx
        # σyy_global = E/(1-ν²) * ν * εxx
        σxx_global = E / (1 - ν^2) * εxx_global
        σyy_global = E / (1 - ν^2) * ν * εxx_global
        
        # Element 1 has local x ≈ global x, local y ≈ global y
        # Element 2 has local x ≈ global y, local y ≈ global -x (rotated 90°)
        # So for element 2: σ_local_xx = σ_global_yy, σ_local_yy = σ_global_xx
        
        # Test element 1 (local ≈ global)
        @test σ1[1] ≈ σxx_global rtol=1e-10  # σxx_local = σxx_global
        @test σ1[2] ≈ σyy_global rtol=1e-10  # σyy_local = σyy_global
        @test σ1[3] ≈ 0.0 atol=1e-6          # τxy = 0
        
        # Test element 2 (local rotated 90° from global)
        @test σ2[1] ≈ σyy_global rtol=1e-10  # σxx_local = σyy_global
        @test σ2[2] ≈ σxx_global rtol=1e-10  # σyy_local = σxx_global
        @test σ2[3] ≈ 0.0 atol=1e-6          # τxy = 0
        
        # Principal stresses should be the same (invariant)
        principal_1 = sort([σ1[1], σ1[2]])
        principal_2 = sort([σ2[1], σ2[2]])
        @test principal_1 ≈ principal_2 rtol=1e-10
        
        @info "Patch test passed" σ1 σ2 σxx_global σyy_global
    end
    
end
