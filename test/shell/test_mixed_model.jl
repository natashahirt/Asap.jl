using Test
using Asap
using LinearAlgebra
using SparseArrays
using Unitful

@testset "Mixed Model Assembly" begin

    @testset "Frame + Shell Assembly" begin
        # Simple structure: two columns with a floor diaphragm
        #
        #   3 ●═══════● 4     5 ●═══════● 6
        #     ║ shell ║         ║ shell ║
        #   1 ●       ● 2       at z=3m (XY plane)
        #     ↑ fixed  
        #
        
        # Create nodes (6 DOF each for frame compatibility)
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)   # fixed base
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :fixed)   # fixed base
        n3 = Asap.Node([0.0u"m", 0.0u"m", 3.0u"m"], :free)    # free top
        n4 = Asap.Node([4.0u"m", 0.0u"m", 3.0u"m"], :free)    # free top
        n5 = Asap.Node([0.0u"m", 4.0u"m", 3.0u"m"], :free)    # floor corner
        n6 = Asap.Node([4.0u"m", 4.0u"m", 3.0u"m"], :free)    # floor corner
        
        nodes = [n1, n2, n3, n4, n5, n6]
        
        # Assign node IDs (6 DOF per node for frame model)
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        n_dof = 36  # 6 nodes × 6 DOF
        
        # Create frame elements (columns)
        sec = Asap.Section(
            0.01u"m^2",   # A
            200e9u"Pa",   # E
            77e9u"Pa",    # G
            1e-4u"m^4",   # Ix
            1e-4u"m^4",   # Iy
            1e-5u"m^4"    # J
        )
        
        col1 = Asap.Element(n1, n3, sec, :column1)
        col2 = Asap.Element(n2, n4, sec, :column2)
        
        # Assign frame element globalIDs
        col1.globalID = [n1.globalID; n3.globalID]
        col2.globalID = [n2.globalID; n4.globalID]
        
        # Process frame elements
        Asap.process_elements!([col1, col2])
        
        # Shell in XY plane at z=3m connecting n3, n4, n6, n5
        shell = Asap.ShellQuad4((n3, n4, n6, n5), 0.15u"m", 30e9u"Pa", 0.2, :floor)
        Asap.process!(shell)
        Asap.populate_globalID!(shell)
        
        # Check shell global IDs are translation DOFs only (indices 1,2,3 of each node)
        @test length(shell.globalID) == 12  # 4 nodes × 3 translation DOFs
        expected_shell_ids = [13, 14, 15, 19, 20, 21, 31, 32, 33, 25, 26, 27]
        @test shell.globalID == expected_shell_ids
        
        # Assemble combined stiffness - convert to common element type
        all_elements = Asap.AbstractElement[col1, col2, shell]
        
        # Test assembly function
        S = Asap.assemble_stiffness(all_elements, n_dof)
        
        @test size(S) == (n_dof, n_dof)
        @test S ≈ S' rtol=1e-10  # symmetric
        
        # Check that shell contributes to translation DOFs
        shell_dofs = shell.globalID
        @test any(S[shell_dofs, shell_dofs] .!= 0)
        
        @info "Mixed model assembly test passed" n_frame_elements=2 n_shell_elements=1 total_dof=n_dof
    end
    
    @testset "Simple Diaphragm Effect" begin
        # Two columns with a floor shell connecting their tops
        # Apply lateral load, check that diaphragm ties them together
        #
        #     P →  3 ●═══════● 4
        #           ║ shell  ║
        #         1 ●        ● 2
        #           ↑ fixed  ↑
        
        # Nodes
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([0.0u"m", 0.0u"m", 3.0u"m"], :free)
        n4 = Asap.Node([4.0u"m", 0.0u"m", 3.0u"m"], :free)
        
        nodes = [n1, n2, n3, n4]
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        n_dof = 24
        
        # Columns
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                          1e-4u"m^4", 1e-4u"m^4", 1e-5u"m^4")
        col1 = Asap.Element(n1, n3, sec)
        col2 = Asap.Element(n2, n4, sec)
        
        # Assign frame element globalIDs
        col1.globalID = [n1.globalID; n3.globalID]
        col2.globalID = [n2.globalID; n4.globalID]
        
        # Process
        Asap.process_elements!([col1, col2])
        
        # Assemble
        elements = Asap.AbstractElement[col1, col2]
        S = Asap.assemble_stiffness(elements, n_dof)
        
        # Apply lateral load at node 3
        F = zeros(n_dof)
        F[13] = 1000.0  # 1 kN in x-direction at node 3
        
        # Fixed DOFs
        fixed_dofs = vcat(collect(n1.globalID), collect(n2.globalID))
        free_dofs = setdiff(1:n_dof, fixed_dofs)
        
        # Solve
        u = zeros(n_dof)
        u[free_dofs] = S[free_dofs, free_dofs] \ F[free_dofs]
        
        # Node 3 should deflect in x, node 4 should not (no coupling without diaphragm)
        u3_x = u[13]
        u4_x = u[19]
        
        @test u3_x > 0  # loaded node deflects
        @test abs(u4_x) < abs(u3_x) * 0.01  # unloaded node doesn't move much
        
        @info "Frame-only baseline" u3_x u4_x
    end
    
    @testset "Rigid Diaphragm Coupling" begin
        # Four columns with a floor diaphragm - load at one corner couples all
        #
        #  P→ 5 ●═══════● 6    (floor at z=3m)
        #       ║       ║
        #     7 ●═══════● 8    (floor at z=3m)
        #
        #     1 ●       ● 2    (base at z=0, fixed)
        #     3 ●       ● 4    (base at z=0, fixed)
        
        # Base nodes (fixed)
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([0.0u"m", 4.0u"m", 0.0u"m"], :fixed)
        n4 = Asap.Node([4.0u"m", 4.0u"m", 0.0u"m"], :fixed)
        
        # Top nodes (free) - at z=3m
        n5 = Asap.Node([0.0u"m", 0.0u"m", 3.0u"m"], :free)
        n6 = Asap.Node([4.0u"m", 0.0u"m", 3.0u"m"], :free)
        n7 = Asap.Node([0.0u"m", 4.0u"m", 3.0u"m"], :free)
        n8 = Asap.Node([4.0u"m", 4.0u"m", 3.0u"m"], :free)
        
        nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
        for (i, node) in enumerate(nodes)
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        n_dof = 48  # 8 nodes × 6 DOF
        
        # Four columns
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                          1e-4u"m^4", 1e-4u"m^4", 1e-5u"m^4")
        col1 = Asap.Element(n1, n5, sec, release=:fixedfixed)
        col2 = Asap.Element(n2, n6, sec, release=:fixedfixed)
        col3 = Asap.Element(n3, n7, sec, release=:fixedfixed)
        col4 = Asap.Element(n4, n8, sec, release=:fixedfixed)
        
        # Assign frame element globalIDs
        col1.globalID = [n1.globalID; n5.globalID]
        col2.globalID = [n2.globalID; n6.globalID]
        col3.globalID = [n3.globalID; n7.globalID]
        col4.globalID = [n4.globalID; n8.globalID]
        
        # Process frames
        Asap.process_elements!([col1, col2, col3, col4])
        
        # Rigid diaphragm shell (very high stiffness) in XY plane at z=3m
        # Quad connecting n5-n6-n8-n7
        rigid_shell = Asap.ShellQuad4((n5, n6, n8, n7), 0.2u"m", 1e12u"Pa", 0.001, :rigid)
        Asap.process!(rigid_shell)
        Asap.populate_globalID!(rigid_shell)
        
        # Assemble combined stiffness
        elements = Asap.AbstractElement[col1, col2, col3, col4, rigid_shell]
        S = Asap.assemble_stiffness(elements, n_dof)
        
        # Apply lateral load at node 5 (corner)
        F = zeros(n_dof)
        F[n5.globalID[1]] = 1000.0  # 1 kN in x-direction
        
        # Fixed DOFs from node definitions
        fixed_dofs = Int[]
        for node in nodes
            for j in 1:6
                if !node.dof[j]
                    push!(fixed_dofs, node.globalID[j])
                end
            end
        end
        free_dofs = setdiff(1:n_dof, fixed_dofs)
        
        # Solve
        u = zeros(n_dof)
        K_ff = Matrix(S[free_dofs, free_dofs])
        u[free_dofs] = K_ff \ F[free_dofs]
        
        # Get x-displacements of all top nodes
        u5_x = u[n5.globalID[1]]
        u6_x = u[n6.globalID[1]]
        u7_x = u[n7.globalID[1]]
        u8_x = u[n8.globalID[1]]
        
        # Diaphragm should couple all top nodes
        @test u5_x > 0  # loaded node deflects
        @test u6_x > 0  # adjacent node deflects
        @test u7_x > 0  # side node deflects
        @test u8_x > 0  # far corner deflects
        
        # Key test: nodes on same row should move together (membrane shell)
        # Front row (5,6) should have nearly equal x-displacement
        @test isapprox(u5_x, u6_x, rtol=0.01)  # same row coupling
        # Back row (7,8) should also have nearly equal x-displacement
        @test isapprox(u7_x, u8_x, rtol=0.01)  # same row coupling
        
        # The coupling ratio shows how much the diaphragm transfers load
        # For stiff diaphragm, back row should deflect > 25% of front row
        coupling_ratio = u7_x / u5_x
        @test coupling_ratio > 0.25  # significant load transfer
        
        @info "Rigid diaphragm coupling" u5_x u6_x u7_x u8_x coupling_ratio
    end

end
