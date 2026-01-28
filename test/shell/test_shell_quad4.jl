using Test
using Asap
using LinearAlgebra
using SparseArrays
using Unitful

@testset "ShellQuad4 Element" begin

    @testset "Geometry & LCS" begin
        # Simple square in XY plane
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        elem = Asap.ShellQuad4((n1, n2, n3, n4), 0.1u"m", 200e9u"Pa", 0.3)
        Asap.lcs!(elem)
        
        # Check LCS
        x_loc, y_loc, z_loc = elem.LCS
        @test x_loc ≈ [1.0, 0.0, 0.0] atol=1e-10
        @test y_loc ≈ [0.0, 1.0, 0.0] atol=1e-10
        @test z_loc ≈ [0.0, 0.0, 1.0] atol=1e-10
        
        # Check area (1×1 square)
        @test elem.area ≈ 1.0 atol=1e-10
    end
    
    @testset "Stiffness Matrix Properties" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        elem = Asap.ShellQuad4((n1, n2, n3, n4), 0.1u"m", 200e9u"Pa", 0.3)
        Asap.process!(elem)
        
        # Local K: 8×8, symmetric
        K_local = Asap.local_K(elem)
        @test size(K_local) == (8, 8)
        @test K_local ≈ K_local' rtol=1e-10
        
        # Global K: 12×12, symmetric
        @test size(elem.K) == (12, 12)
        @test elem.K ≈ elem.K' rtol=1e-10
        
        # Eigenvalues: 8 DOFs - 3 RBM = 5 non-zero for local
        # Global 12×12 has rank ≤ 5, so 7+ zero eigenvalues
        eigvals_K = eigvals(Symmetric(elem.K))
        @test all(eigvals_K .>= -1e-6 * maximum(eigvals_K))
    end
    
    @testset "Patch Test - Constant Strain" begin
        # Single quad element, apply uniform strain εxx
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        nodes = [n1, n2, n3, n4]
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*3+1 : i*3)
        end
        
        elem = Asap.ShellQuad4((n1, n2, n3, n4), 0.1u"m", 200e9u"Pa", 0.3)
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in nodes]...)
        
        # Apply uniform x-strain
        εxx = 0.001
        u = zeros(12)
        u[1] = 0.0        # n1 x
        u[4] = εxx * 1.0  # n2 x
        u[7] = εxx * 1.0  # n3 x
        u[10] = 0.0       # n4 x
        
        σ = Asap.stress(elem, u)
        
        E = 200e9
        ν = 0.3
        σxx_expected = E / (1 - ν^2) * εxx
        σyy_expected = E / (1 - ν^2) * ν * εxx
        
        @test σ[1] ≈ σxx_expected rtol=1e-10
        @test σ[2] ≈ σyy_expected rtol=1e-10
        @test σ[3] ≈ 0.0 atol=1e-6
    end
    
    @testset "Cook's Membrane Benchmark" begin
        # Cook's tapered membrane: classic FEA benchmark
        # Reference: Cook, Malkus, Plesha "Concepts and Applications of FEA"
        # 
        #            ←─ P (distributed shear)
        #         ●───●
        #        /    │   44 units height at left
        #       /     │
        #      /      │
        #     /       │  16 units height at right
        #    ●────────●
        #        48 units width
        #    ↑ fixed edge (left)
        #
        # Standard parameters (dimensionless/normalized):
        #   E = 1 (unit modulus, effectively MPa = N/mm²)
        #   ν = 1/3
        #   t = 1 (unit thickness)
        #   P = 1 N total shear on right edge
        # Reference solution: tip vertical displacement ≈ 23.9 (with fine mesh)
        # Single Q4 element typically gives ~11-15 due to locking in distorted quad
        
        # Node positions (using meters internally but scaled appropriately)
        # We use E=1e6 Pa (1 MPa) with mm dimensions for consistency
        n1 = Asap.Node([0.0u"mm", 0.0u"mm", 0.0u"mm"], :free)
        n2 = Asap.Node([48.0u"mm", 0.0u"mm", 0.0u"mm"], :free)
        n3 = Asap.Node([48.0u"mm", 16.0u"mm", 0.0u"mm"], :free)
        n4 = Asap.Node([0.0u"mm", 44.0u"mm", 0.0u"mm"], :free)
        
        nodes = [n1, n2, n3, n4]
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*3+1 : i*3)
        end
        
        # Standard benchmark parameters
        E = 1.0e6  # 1 MPa = 1 N/mm² (unit modulus in mm system)
        ν = 1/3
        t = 1.0    # 1 mm thickness
        
        elem = Asap.ShellQuad4((n1, n2, n3, n4), t*u"mm", E*u"Pa", ν)
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in nodes]...)
        
        # Assemble (single element)
        n_dof = 12
        S = zeros(n_dof, n_dof)
        idx = elem.globalID
        for i in 1:12, j in 1:12
            S[idx[i], idx[j]] += elem.K[i, j]
        end
        
        # Boundary conditions: 
        # - Left edge (nodes 1, 4) fixed: x, y (and z for 2D constraint)
        # - All z DOFs fixed (2D in-plane problem)
        # Load: distributed shear P=1 N on right edge (nodes 2, 3)
        P_total = 1.0  # N (total force)
        F = zeros(n_dof)
        F[5] = P_total / 2   # n2 y-direction
        F[8] = P_total / 2   # n3 y-direction
        
        # Fixed DOFs: 
        # nodes 1, 4: x,y,z fixed (DOFs 1,2,3 and 10,11,12)
        # all z DOFs (3, 6, 9, 12) for 2D constraint
        fixed = unique([1, 2, 3, 10, 11, 12, 6, 9])
        free = setdiff(1:12, fixed)
        
        # Solve
        u = zeros(n_dof)
        u[free] = S[free, free] \ F[free]
        
        # Tip displacement (average of nodes 2 and 3 y-displacement)
        # Result is in meters (Asap internal units), convert to mm
        tip_disp_y_m = (u[5] + u[8]) / 2
        tip_disp_y_mm = tip_disp_y_m * 1000  # convert m to mm
        
        # For a single distorted quad element:
        # - Reference (fine mesh): ~23.9 mm
        # - Single Q4 with shear locking: typically ~11-15 mm
        @test tip_disp_y_mm > 5.0   # Should be positive and significant
        @test tip_disp_y_mm < 30.0  # Shouldn't exceed fine mesh reference
        
        @info "Cook's membrane (1 element)" tip_displacement_mm=tip_disp_y_mm reference_fine_mesh=23.9
    end
    
    @testset "Cook's Membrane - Mesh Convergence" begin
        # Test mesh convergence for Cook's membrane benchmark
        # Reference solution: ~23.9 mm tip displacement (enhanced/higher-order elements)
        # 
        # Note: Standard Q4 membrane elements exhibit "shear locking" in bending-dominated
        # problems, causing slower convergence. The reference 23.9 is typically achieved with:
        # - Incompatible modes / enhanced strain
        # - Reduced integration with hourglass control
        # - Higher-order elements (Q8, Q9)
        #
        # Our standard Q4 will converge but requires finer meshes to approach reference.
        
        """
        Create NxN mesh of Cook's membrane and return tip y-displacement in mm.
        
        Geometry: tapered membrane
        - Left edge: x=0, y ∈ [0, 44] mm
        - Right edge: x=48, y ∈ [0, 16] mm
        - Linear interpolation of height across width
        """
        function solve_cook_mesh(N::Int)
            # Standard parameters
            E = 1.0e6  # 1 MPa
            ν = 1/3
            t = 1.0    # mm
            P_total = 1.0  # N total shear
            
            # Geometry bounds
            width = 48.0   # mm
            h_left = 44.0  # mm (left edge height)
            h_right = 16.0 # mm (right edge height)
            
            # Generate node grid: (N+1) x (N+1)
            nodes = Asap.Node[]
            node_map = zeros(Int, N+1, N+1)  # (i,j) → node index
            
            for i in 0:N
                for j in 0:N
                    # Parametric coordinates
                    ξ = i / N  # 0 to 1 across width
                    η = j / N  # 0 to 1 up height
                    
                    # Physical coordinates (linear interpolation)
                    x = ξ * width
                    h_at_x = h_left + ξ * (h_right - h_left)  # height at this x
                    y = η * h_at_x
                    
                    node = Asap.Node([x*u"mm", y*u"mm", 0.0u"mm"], :free)
                    push!(nodes, node)
                    node_map[i+1, j+1] = length(nodes)
                end
            end
            
            # Assign global DOF indices
            for (idx, node) in enumerate(nodes)
                node.nodeID = idx
                node.globalID = collect((idx-1)*3+1 : idx*3)
            end
            
            # Create quad elements
            elements = Asap.ShellQuad4[]
            for i in 1:N
                for j in 1:N
                    # Nodes of this quad (CCW order)
                    n1_idx = node_map[i, j]
                    n2_idx = node_map[i+1, j]
                    n3_idx = node_map[i+1, j+1]
                    n4_idx = node_map[i, j+1]
                    
                    elem = Asap.ShellQuad4(
                        (nodes[n1_idx], nodes[n2_idx], nodes[n3_idx], nodes[n4_idx]),
                        t*u"mm", E*u"Pa", ν
                    )
                    Asap.process!(elem)
                    elem.globalID = vcat([nodes[n1_idx].globalID, nodes[n2_idx].globalID,
                                          nodes[n3_idx].globalID, nodes[n4_idx].globalID]...)
                    push!(elements, elem)
                end
            end
            
            # Assemble global stiffness matrix
            n_nodes = length(nodes)
            n_dof = 3 * n_nodes
            S = zeros(n_dof, n_dof)
            
            for elem in elements
                idx = elem.globalID
                for ii in 1:12, jj in 1:12
                    S[idx[ii], idx[jj]] += elem.K[ii, jj]
                end
            end
            
            # Load: distributed shear on right edge (i = N+1)
            F = zeros(n_dof)
            right_edge_nodes = [node_map[N+1, j] for j in 1:N+1]
            force_per_node = P_total / length(right_edge_nodes)
            for node_idx in right_edge_nodes
                y_dof = (node_idx - 1) * 3 + 2  # y-direction DOF
                F[y_dof] = force_per_node
            end
            
            # Boundary conditions:
            # - Left edge (i=1) x,y fixed
            # - All z DOFs fixed (2D problem)
            fixed_dofs = Int[]
            for j in 1:N+1
                left_node = node_map[1, j]
                push!(fixed_dofs, (left_node-1)*3 + 1)  # x
                push!(fixed_dofs, (left_node-1)*3 + 2)  # y
            end
            for node_idx in 1:n_nodes
                push!(fixed_dofs, (node_idx-1)*3 + 3)  # z
            end
            fixed_dofs = unique(fixed_dofs)
            free_dofs = setdiff(1:n_dof, fixed_dofs)
            
            # Solve
            u = zeros(n_dof)
            u[free_dofs] = S[free_dofs, free_dofs] \ F[free_dofs]
            
            # Tip displacement: average y-displacement of right edge top node
            tip_node = node_map[N+1, N+1]  # top-right corner
            tip_y_dof = (tip_node - 1) * 3 + 2
            tip_disp_mm = u[tip_y_dof] * 1000  # m to mm
            
            return tip_disp_mm
        end
        
        # Test convergence with increasing mesh refinement
        results = Dict{Int, Float64}()
        for N in [1, 2, 4, 8]
            results[N] = solve_cook_mesh(N)
        end
        
        # Verify monotonic convergence (key test for correctness)
        @test results[2] > results[1]  # 2x2 better than 1x1
        @test results[4] > results[2]  # 4x4 better than 2x2
        @test results[8] > results[4]  # 8x8 better than 4x4
        
        # Standard Q4 with shear locking converges slowly
        # 8x8 mesh typically achieves ~65-70% of reference (15-17 mm)
        @test results[8] > 14.0   # Should exceed 14 mm
        @test results[8] < 20.0   # Standard Q4 won't reach 20 mm with 8x8
        
        # Reference is ~23.9 mm (with enhanced elements or very fine mesh)
        @info "Cook's membrane convergence (standard Q4 with shear locking)" mesh_1x1=results[1] mesh_2x2=results[2] mesh_4x4=results[4] mesh_8x8=results[8] reference=23.9
    end
    
    @testset "Cook's Membrane - Triangles vs Quads vs Fine Mesh" begin
        # Compare element types and mesh density for Cook's membrane
        # This tests our Tri3 vs Quad4 and explores convergence behavior
        
        """Solve Cook's membrane with NxN quads."""
        function solve_cook_quads(N::Int)
            E, ν, t, P_total = 1.0e6, 1/3, 1.0, 1.0
            width, h_left, h_right = 48.0, 44.0, 16.0
            
            nodes = Asap.Node[]
            node_map = zeros(Int, N+1, N+1)
            
            for i in 0:N, j in 0:N
                ξ, η = i/N, j/N
                x = ξ * width
                y = η * (h_left + ξ * (h_right - h_left))
                push!(nodes, Asap.Node([x*u"mm", y*u"mm", 0.0u"mm"], :free))
                node_map[i+1, j+1] = length(nodes)
            end
            
            for (idx, node) in enumerate(nodes)
                node.nodeID = idx
                node.globalID = collect((idx-1)*3+1 : idx*3)
            end
            
            elements = Asap.ShellQuad4[]
            for i in 1:N, j in 1:N
                n1, n2 = node_map[i, j], node_map[i+1, j]
                n3, n4 = node_map[i+1, j+1], node_map[i, j+1]
                elem = Asap.ShellQuad4((nodes[n1], nodes[n2], nodes[n3], nodes[n4]), t*u"mm", E*u"Pa", ν)
                Asap.process!(elem)
                elem.globalID = vcat([nodes[n1].globalID, nodes[n2].globalID, nodes[n3].globalID, nodes[n4].globalID]...)
                push!(elements, elem)
            end
            
            n_dof = 3 * length(nodes)
            S, F = zeros(n_dof, n_dof), zeros(n_dof)
            
            for elem in elements
                idx = elem.globalID
                for ii in 1:12, jj in 1:12
                    S[idx[ii], idx[jj]] += elem.K[ii, jj]
                end
            end
            
            for j in 1:N+1
                F[(node_map[N+1, j]-1)*3 + 2] = P_total / (N+1)
            end
            
            fixed = unique(vcat(
                [(node_map[1,j]-1)*3+k for j in 1:N+1 for k in 1:2],
                [(i-1)*3+3 for i in 1:length(nodes)]
            ))
            free = setdiff(1:n_dof, fixed)
            
            u = zeros(n_dof)
            u[free] = S[free, free] \ F[free]
            return u[(node_map[N+1, N+1]-1)*3 + 2] * 1000
        end
        
        """Solve Cook's membrane with NxN grid split into 2*N*N triangles."""
        function solve_cook_triangles(N::Int)
            E, ν, t, P_total = 1.0e6, 1/3, 1.0, 1.0
            width, h_left, h_right = 48.0, 44.0, 16.0
            
            nodes = Asap.Node[]
            node_map = zeros(Int, N+1, N+1)
            
            for i in 0:N, j in 0:N
                ξ, η = i/N, j/N
                x = ξ * width
                y = η * (h_left + ξ * (h_right - h_left))
                push!(nodes, Asap.Node([x*u"mm", y*u"mm", 0.0u"mm"], :free))
                node_map[i+1, j+1] = length(nodes)
            end
            
            for (idx, node) in enumerate(nodes)
                node.nodeID = idx
                node.globalID = collect((idx-1)*3+1 : idx*3)
            end
            
            # Each quad cell split into 2 triangles (diagonal split)
            elements = Asap.ShellTri3[]
            for i in 1:N, j in 1:N
                n1, n2 = node_map[i, j], node_map[i+1, j]
                n3, n4 = node_map[i+1, j+1], node_map[i, j+1]
                
                # Lower triangle: n1-n2-n3
                elem1 = Asap.ShellTri3((nodes[n1], nodes[n2], nodes[n3]), t*u"mm", E*u"Pa", ν)
                Asap.process!(elem1)
                elem1.globalID = vcat([nodes[n1].globalID, nodes[n2].globalID, nodes[n3].globalID]...)
                push!(elements, elem1)
                
                # Upper triangle: n1-n3-n4
                elem2 = Asap.ShellTri3((nodes[n1], nodes[n3], nodes[n4]), t*u"mm", E*u"Pa", ν)
                Asap.process!(elem2)
                elem2.globalID = vcat([nodes[n1].globalID, nodes[n3].globalID, nodes[n4].globalID]...)
                push!(elements, elem2)
            end
            
            n_dof = 3 * length(nodes)
            S, F = zeros(n_dof, n_dof), zeros(n_dof)
            
            for elem in elements
                idx = elem.globalID
                for ii in 1:9, jj in 1:9
                    S[idx[ii], idx[jj]] += elem.K[ii, jj]
                end
            end
            
            for j in 1:N+1
                F[(node_map[N+1, j]-1)*3 + 2] = P_total / (N+1)
            end
            
            fixed = unique(vcat(
                [(node_map[1,j]-1)*3+k for j in 1:N+1 for k in 1:2],
                [(i-1)*3+3 for i in 1:length(nodes)]
            ))
            free = setdiff(1:n_dof, fixed)
            
            u = zeros(n_dof)
            u[free] = S[free, free] \ F[free]
            return u[(node_map[N+1, N+1]-1)*3 + 2] * 1000
        end
        
        # Compare at different resolutions
        # Note: 8x8 quads = 64 elements, 8x8 triangles = 128 elements
        quad_8 = solve_cook_quads(8)
        quad_16 = solve_cook_quads(16)
        tri_8 = solve_cook_triangles(8)
        tri_16 = solve_cook_triangles(16)
        
        # Triangles (CST) are even stiffer than Q4 due to constant strain assumption
        # So we expect: tri < quad at same node count
        @test tri_8 < quad_8   # Tri3 stiffer than Q4 at same grid
        @test tri_16 < quad_16
        
        # Both should converge with refinement
        @test quad_16 > quad_8
        @test tri_16 > tri_8
        
        # Q4 should be closer to reference than Tri3 at same node density
        @test quad_16 > tri_16
        
        @info "Cook's membrane: Triangles vs Quads" quad_8x8=quad_8 quad_16x16=quad_16 tri_8x8=tri_8 tri_16x16=tri_16 reference=23.9
    end
    
end
