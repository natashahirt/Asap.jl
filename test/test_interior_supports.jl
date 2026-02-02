#=
Interior Supports Tests
=======================

Tests for Shell() with interior_supports parameter.
Verifies correct meshing, fixity application, and moment distribution.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "Interior Supports" begin
    
    @testset "Basic interior support creation" begin
        # 10m x 6m slab with interior beam at y=3m
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([10.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([10.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        
        # Interior beam from (0,3) to (10,3)
        beam_n1 = Node([0.0u"m", 3.0u"m", 0.0u"m"], :pinned)
        beam_n2 = Node([10.0u"m", 3.0u"m", 0.0u"m"], :pinned)
        beam_sec = Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        interior_beam = Element(beam_n1, beam_n2, beam_sec)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        # Create shells with interior support
        shells = Shell((n1, n2, n3, n4), sec;
            interior_supports = [interior_beam],
            edge_support_type = :pinned,
            interior_support_type = :pinned
        )
        
        @test length(shells) >= 4  # Should have multiple shell elements
        
        # Check that interior support nodes exist
        all_nodes = get_nodes(shells)
        interior_support_nodes = filter(n -> begin
            y = ustrip(u"m", n.position[2])
            abs(y - 3.0) < 0.01
        end, all_nodes)
        
        @test length(interior_support_nodes) >= 3  # Should have nodes along y=3m
    end
    
    @testset "Node fixity application" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        # Interior beam at y=2m
        beam_n1 = Node([0.0u"m", 2.0u"m", 0.0u"m"], :pinned)
        beam_n2 = Node([4.0u"m", 2.0u"m", 0.0u"m"], :pinned)
        beam_sec = Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        interior_beam = Element(beam_n1, beam_n2, beam_sec)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        shells = Shell((n1, n2, n3, n4), 4, sec;
            interior_supports = [interior_beam],
            edge_support_type = :pinned,
            interior_support_type = :zfixed  # Different fixity for interior
        )
        
        all_nodes = get_nodes(shells)
        
        # Check edge nodes (not corners) are :pinned
        edge_nodes = filter(n -> begin
            x, y = ustrip(u"m", n.position[1]), ustrip(u"m", n.position[2])
            on_edge = (abs(x) < 0.01 || abs(x - 4) < 0.01 || abs(y) < 0.01 || abs(y - 4) < 0.01)
            not_corner = !((abs(x) < 0.01 || abs(x - 4) < 0.01) && (abs(y) < 0.01 || abs(y - 4) < 0.01))
            not_interior_support = abs(y - 2.0) > 0.01
            on_edge && not_corner && not_interior_support
        end, all_nodes)
        
        # Edge nodes should have pinned fixity (xyz fixed, rotations free)
        for node in edge_nodes
            @test node.dof[1:3] == [false, false, false]  # xyz fixed
            @test node.dof[4:6] == [true, true, true]     # rotations free
        end
        
        # Check interior support nodes have :zfixed
        interior_support_nodes = filter(n -> begin
            y = ustrip(u"m", n.position[2])
            x = ustrip(u"m", n.position[1])
            abs(y - 2.0) < 0.01 && x > 0.01 && x < 3.99  # Interior support, not corners
        end, all_nodes)
        
        for node in interior_support_nodes
            @test node.dof[3] == false  # z fixed
            @test node.dof[1] == true   # x free
            @test node.dof[2] == true   # y free
        end
    end
    
    @testset "Node pair support" begin
        # Test using (Node, Node) pairs instead of Elements
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        support_start = Node([0.0u"m", 2.0u"m", 0.0u"m"], :pinned)
        support_end = Node([4.0u"m", 2.0u"m", 0.0u"m"], :pinned)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        shells = Shell((n1, n2, n3, n4), sec;
            interior_supports = [(support_start, support_end)]
        )
        
        all_nodes = get_nodes(shells)
        
        # Should have nodes along y=2m
        interior_nodes = filter(n -> abs(ustrip(u"m", n.position[2]) - 2.0) < 0.01, all_nodes)
        @test length(interior_nodes) >= 3
    end
    
    @testset "Moment distribution - single span vs two-span" begin
        # Compare deflection of single span vs two-span continuous slab
        # Two-span should have LESS max deflection due to continuity
        
        L = 6.0  # meters
        w = 4.0  # meters
        pressure = 5.0u"kPa"
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
        
        # --- Single span slab (no interior support) ---
        n1_ss = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2_ss = Node([L*u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3_ss = Node([L*u"m", w*u"m", 0.0u"m"], :pinned)
        n4_ss = Node([0.0u"m", w*u"m", 0.0u"m"], :pinned)
        
        shells_ss = Shell((n1_ss, n2_ss, n3_ss, n4_ss), 6, sec;
            edge_support_type = :pinned
        )
        
        loads_ss = [AreaLoad(shells_ss, pressure)]
        model_ss = ShellModel(get_nodes(shells_ss), shells_ss, loads_ss)
        solve!(model_ss)
        
        # Find max deflection
        all_nodes_ss = get_nodes(shells_ss)
        max_defl_ss = maximum(abs(ustrip(u"m", n.displacement[3])) for n in all_nodes_ss)
        
        # --- Two-span slab (interior support at midspan) ---
        n1_ts = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2_ts = Node([L*u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3_ts = Node([L*u"m", w*u"m", 0.0u"m"], :pinned)
        n4_ts = Node([0.0u"m", w*u"m", 0.0u"m"], :pinned)
        
        # Interior support at x = L/2
        support_n1 = Node([(L/2)*u"m", 0.0u"m", 0.0u"m"], :pinned)
        support_n2 = Node([(L/2)*u"m", w*u"m", 0.0u"m"], :pinned)
        support_sec = Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        interior_beam = Element(support_n1, support_n2, support_sec)
        
        shells_ts = Shell((n1_ts, n2_ts, n3_ts, n4_ts), 6, sec;
            interior_supports = [interior_beam],
            edge_support_type = :pinned,
            interior_support_type = :pinned
        )
        
        loads_ts = [AreaLoad(shells_ts, pressure)]
        model_ts = ShellModel(get_nodes(shells_ts), shells_ts, loads_ts)
        solve!(model_ts)
        
        all_nodes_ts = get_nodes(shells_ts)
        max_defl_ts = maximum(abs(ustrip(u"m", n.displacement[3])) for n in all_nodes_ts)
        
        # Interior support nodes should have zero deflection
        interior_support_nodes = filter(n -> begin
            x = ustrip(u"m", n.position[1])
            abs(x - L/2) < 0.01
        end, all_nodes_ts)
        
        for node in interior_support_nodes
            @test abs(ustrip(u"m", node.displacement[3])) < 1e-10
        end
        
        # Two-span should have MUCH less max deflection than single span
        # For simply-supported: ratio should be roughly (1/2)^4 = 1/16 for uniform load
        # But with continuity effects, it's more like 1/5 to 1/8
        @test max_defl_ts < max_defl_ss / 3  # Conservative check
        
        @info "Deflection comparison" single_span=max_defl_ss two_span=max_defl_ts ratio=max_defl_ss/max_defl_ts
    end
    
    @testset "Free edges" begin
        # Test that edge_support_type = :free leaves edges unconstrained
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)  # Only corners fixed
        n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Node([2.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        n4 = Node([0.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        shells = Shell((n1, n2, n3, n4), 4, sec;
            edge_support_type = :free  # Edges not supported
        )
        
        all_nodes = get_nodes(shells)
        
        # Non-corner edge nodes should be free
        edge_nodes = filter(n -> begin
            x, y = ustrip(u"m", n.position[1]), ustrip(u"m", n.position[2])
            on_edge = (abs(x) < 0.01 || abs(x - 2) < 0.01 || abs(y) < 0.01 || abs(y - 2) < 0.01)
            not_corner = !((abs(x) < 0.01 || abs(x - 2) < 0.01) && (abs(y) < 0.01 || abs(y - 2) < 0.01))
            on_edge && not_corner
        end, all_nodes)
        
        for node in edge_nodes
            @test all(node.dof)  # All DOFs should be free
        end
    end
    
    @testset "Multiple interior supports" begin
        # Test multiple interior beams
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([9.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([9.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        # Two interior beams at x=3m and x=6m
        beam_sec = Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        beam1_n1 = Node([3.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        beam1_n2 = Node([3.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        beam1 = Element(beam1_n1, beam1_n2, beam_sec)
        
        beam2_n1 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        beam2_n2 = Node([6.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        beam2 = Element(beam2_n1, beam2_n2, beam_sec)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        
        shells = Shell((n1, n2, n3, n4), 6, sec;
            interior_supports = [beam1, beam2]
        )
        
        all_nodes = get_nodes(shells)
        
        # Should have nodes along both support lines
        support1_nodes = filter(n -> abs(ustrip(u"m", n.position[1]) - 3.0) < 0.01, all_nodes)
        support2_nodes = filter(n -> abs(ustrip(u"m", n.position[1]) - 6.0) < 0.01, all_nodes)
        
        @test length(support1_nodes) >= 3
        @test length(support2_nodes) >= 3
    end
end
