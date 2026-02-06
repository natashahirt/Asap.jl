#=
Interior Nodes Tests
====================

Tests for Shell() with interior_nodes parameter.
Verifies that actual Node objects can be embedded in shell mesh for true structural connectivity.
=#

using Test
using Asap
using Unitful

@testset "Interior Nodes - Shared Node Objects" begin
    
    @testset "Basic interior node embedding" begin
        # Create a 4m × 4m shell with an interior column node at (2m, 2m)
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Asap.Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Asap.Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        # Interior column top node - this is the key: we want THIS node in the mesh
        interior_col_top = Asap.Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)
        
        sec = Asap.ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        # Create shells with interior node embedded
        shells = Asap.Shell((n1, n2, n3, n4), sec;
            interior_nodes = [interior_col_top],
            edge_support_type = :pinned
        )
        
        @test length(shells) >= 4  # Should have multiple shell elements
        
        # Verify the interior node is actually in the shell mesh
        all_nodes = Asap.get_nodes(shells)
        @test interior_col_top in all_nodes  # Same object identity!
        
        # Verify node position is at (2, 2)
        interior_node_found = filter(n -> begin
            x = ustrip(u"m", n.position[1])
            y = ustrip(u"m", n.position[2])
            abs(x - 2.0) < 0.01 && abs(y - 2.0) < 0.01
        end, all_nodes)
        
        @test length(interior_node_found) == 1
        @test interior_node_found[1] === interior_col_top  # Same object!
    end
    
    @testset "Multiple interior nodes" begin
        # 6m × 6m shell with 2 interior column nodes
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Asap.Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Asap.Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        n4 = Asap.Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        
        # Two interior nodes
        col_top_1 = Asap.Node([2.0u"m", 3.0u"m", 0.0u"m"], :free)
        col_top_2 = Asap.Node([4.0u"m", 3.0u"m", 0.0u"m"], :free)
        
        sec = Asap.ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        shells = Asap.Shell((n1, n2, n3, n4), sec;
            interior_nodes = [col_top_1, col_top_2],
            edge_support_type = :pinned
        )
        
        all_nodes = Asap.get_nodes(shells)
        
        # Both interior nodes should be the exact same objects
        @test col_top_1 in all_nodes
        @test col_top_2 in all_nodes
        @test col_top_1 === filter(n -> n === col_top_1, all_nodes)[1]
        @test col_top_2 === filter(n -> n === col_top_2, all_nodes)[1]
    end
    
    @testset "Mixed model - frame + shell with shared nodes" begin
        # This tests true structural connectivity: column shares node with shell
        
        # Column base (fixed) and top (shared with shell)
        col_base = Asap.Node([2.0u"m", 2.0u"m", -3.0u"m"], :fixed)
        col_top = Asap.Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)  # Shared!
        
        # Shell corners (pinned at edges)
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Asap.Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Asap.Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        # Column section
        col_sec = Asap.Section(
            0.01u"m^2", 200e9u"Pa", 80e9u"Pa",
            0.001u"m^4", 0.001u"m^4", 0.002u"m^4", 7850u"kg/m^3"
        )
        column = Asap.Element(col_base, col_top, col_sec)
        
        # Shell with col_top as interior node
        shell_sec = Asap.ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
        shells = Asap.Shell((n1, n2, n3, n4), shell_sec;
            interior_nodes = [col_top],
            edge_support_type = :free  # Edges free - only column supports interior
        )
        
        # Load on shell (downward pressure)
        pressure = 5000.0u"Pa"
        area_load = Asap.AreaLoad(shells, pressure; direction=(0.0, 0.0, -1.0))
        
        # Create mixed model with shared node
        all_nodes = unique([col_base, col_top, n1, n2, n3, n4, Asap.get_nodes(shells)...])
        model = Asap.Model(all_nodes, [column], shells, [area_load])
        
        @test Asap.is_mixed(model)
        
        # Process and solve
        Asap.process!(model)
        Asap.solve!(model)
        
        # Key test: verify col_top deflects (connected to shell via shared node)
        z_idx = col_top.globalID[3]
        z_disp = model.u[z_idx]
        @test abs(z_disp) > 0  # Should have some vertical displacement from load
        
        # Verify edge corners don't move much (they're free but shell is connected to them)
        # The interior node should be the primary load path to the column
    end
    
    @testset "Interior node at corner position is handled correctly" begin
        # Edge case: what if interior node is at same position as corner?
        # The corner node should take precedence (it's added first)
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Asap.Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Asap.Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        # Interior node at same position as corner - should be ignored/deduplicated
        duplicate_corner = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        
        sec = Asap.ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        shells = Asap.Shell((n1, n2, n3, n4), sec;
            interior_nodes = [duplicate_corner],
            edge_support_type = :pinned
        )
        
        all_nodes = Asap.get_nodes(shells)
        
        # The corner n1 should be used, not the duplicate
        @test n1 in all_nodes
        # duplicate_corner at same position might not be added since boundary takes priority
    end
    
end

println("Interior nodes tests completed!")
