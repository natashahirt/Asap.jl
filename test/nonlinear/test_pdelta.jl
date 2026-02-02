#=
P-Delta Analysis Tests
======================

Validation cases for P-delta iteration.
=#

using Test
using Asap
using Unitful

@testset "P-Delta Analysis" begin
    
    @testset "Simple Frame P-Delta" begin
        # Two-story frame under gravity + lateral load
        # P-delta should amplify lateral displacement
        
        L_bay = 6.0u"m"
        H_story = 3.5u"m"
        
        # Steel W-section properties (approximate)
        E = 200e9u"Pa"
        G = 77e9u"Pa"
        
        # Column section
        col_section = Section(
            0.015u"m^2", E, G,
            3e-4u"m^4", 1e-4u"m^4", 5e-6u"m^4"
        )
        
        # Beam section (larger)
        beam_section = Section(
            0.02u"m^2", E, G,
            5e-4u"m^4", 2e-4u"m^4", 8e-6u"m^4"
        )
        
        # Nodes: 2 stories, single bay
        nodes = [
            # Ground level (fixed)
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([L_bay, 0.0u"m", 0.0u"m"], :fixed),
            # Level 1
            Node([0.0u"m", 0.0u"m", H_story], :free),
            Node([L_bay, 0.0u"m", H_story], :free),
            # Level 2 (roof)
            Node([0.0u"m", 0.0u"m", 2*H_story], :free),
            Node([L_bay, 0.0u"m", 2*H_story], :free),
        ]
        
        elements = [
            # Columns
            Element(nodes[1], nodes[3], col_section, :col_L1_left),
            Element(nodes[2], nodes[4], col_section, :col_L1_right),
            Element(nodes[3], nodes[5], col_section, :col_L2_left),
            Element(nodes[4], nodes[6], col_section, :col_L2_right),
            # Beams
            Element(nodes[3], nodes[4], beam_section, :beam_L1),
            Element(nodes[5], nodes[6], beam_section, :beam_L2),
        ]
        
        # Gravity + lateral loads
        P_gravity = 500e3u"N"  # 500 kN per column
        V_lateral = 50e3u"N"   # 50 kN lateral at each floor
        
        loads = [
            # Gravity at column tops
            NodeForce(nodes[3], [0.0u"N", 0.0u"N", -P_gravity]),
            NodeForce(nodes[4], [0.0u"N", 0.0u"N", -P_gravity]),
            NodeForce(nodes[5], [0.0u"N", 0.0u"N", -P_gravity]),
            NodeForce(nodes[6], [0.0u"N", 0.0u"N", -P_gravity]),
            # Lateral at floors
            NodeForce(nodes[3], [V_lateral, 0.0u"N", 0.0u"N"]),
            NodeForce(nodes[4], [V_lateral, 0.0u"N", 0.0u"N"]),
            NodeForce(nodes[5], [V_lateral, 0.0u"N", 0.0u"N"]),
            NodeForce(nodes[6], [V_lateral, 0.0u"N", 0.0u"N"]),
        ]
        
        model = Model(nodes, elements, loads)
        
        # First-order analysis
        solve!(model)
        u_first_order = maximum(abs.(model.u))
        
        # P-delta analysis
        result = solve_pdelta!(model; max_iter=10, tol=1e-4, verbose=false)
        u_pdelta = maximum(abs.(model.u))
        
        @test result.converged
        @test result.iterations <= 10
        @test result.amplification > 1.0  # P-delta should amplify
        @test result.amplification < 2.0  # But not by too much for stable frame
        
        # Verify displacement increased
        @test u_pdelta > u_first_order
        
        println("P-delta frame test:")
        println("  First-order max disp: $(round(u_first_order * 1000, digits=2)) mm")
        println("  P-delta max disp: $(round(u_pdelta * 1000, digits=2)) mm")
        println("  Amplification: $(round(result.amplification, digits=3))")
        println("  Iterations: $(result.iterations)")
    end
    
    @testset "P-delta Amplification Behavior" begin
        # Test that P-delta gives reasonable amplification for a frame
        # under combined gravity and lateral load
        
        L = 4.0u"m"
        E = 200e9u"Pa"
        G = 77e9u"Pa"
        
        # Moderately slender section
        section = Section(
            0.01u"m^2", E, G,
            1e-4u"m^4", 1e-4u"m^4", 1e-5u"m^4"
        )
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([0.0u"m", 0.0u"m", L/2], :free),
            Node([0.0u"m", 0.0u"m", L], :free)
        ]
        
        elements = [
            Element(nodes[1], nodes[2], section),
            Element(nodes[2], nodes[3], section)
        ]
        
        # Moderate axial load + lateral load
        P = 200e3u"N"  # Compression
        V = 5e3u"N"    # Lateral
        
        loads = [
            NodeForce(nodes[3], [V, 0.0u"N", -P]),
        ]
        
        model = Model(nodes, elements, loads)
        result = solve_pdelta!(model; max_iter=10, tol=1e-3)
        
        # Should converge with amplification > 1 (P-delta effect)
        @test result.converged
        @test result.amplification >= 1.0
        
        println("Column P-delta:")
        println("  Converged: $(result.converged)")
        println("  Amplification: $(round(result.amplification, digits=3))")
    end
    
    @testset "B2 Factor Convenience" begin
        # Simple test that B2_factor function works
        L = 3.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([0.0u"m", 0.0u"m", L], :free)
        ]
        elements = [Element(nodes[1], nodes[2], section)]
        
        loads = [
            NodeForce(nodes[2], [1e3u"N", 0.0u"N", -100e3u"N"])  # Lateral + gravity
        ]
        
        model = Model(nodes, elements, loads)
        
        B2 = B2_factor(model)
        @test B2 > 1.0
        @test B2 < 10.0  # Reasonable range
    end
end
