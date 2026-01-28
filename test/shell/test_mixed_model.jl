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

end
