#=
Shell-Beam Interaction Tests
============================

Tests for converting shell panel geometry into TributaryLoads for edge beams.

Test Categories:
1. Edge Matching - Foundation for everything
2. Shell Boundary Extraction - Extract boundary from triangle mesh
3. Shell-to-Tributary Conversion - Core function
4. One-Way vs Isotropic Distribution - Verify zero load on non-bearing beams
5. Interior Beam Handling - Error on unmatched beams
6. Multi-Panel Support - Accumulation on shared interior beams
7. Equilibrium Tests - Sanity checks
=#

using Test
using Asap
using Meshes
using Unitful
using LinearAlgebra

# =============================================================================
# Test Utilities
# =============================================================================

"""Compute beam length from node positions (for unprocessed elements)."""
function beam_length_from_nodes(elem)
    dx = ustrip(u"m", elem.nodeEnd.position[1] - elem.nodeStart.position[1])
    dy = ustrip(u"m", elem.nodeEnd.position[2] - elem.nodeStart.position[2])
    dz = ustrip(u"m", elem.nodeEnd.position[3] - elem.nodeStart.position[3])
    return sqrt(dx^2 + dy^2 + dz^2)
end

"""Create a simple rectangular set of shell elements (4 triangles meeting at center)."""
function make_rect_shells(Lx, Ly; z=0.0, thickness=0.15, E=30e9, ν=0.2, ρ=2400.0)
    # Corner nodes
    n1 = Asap.Node([0.0u"m", 0.0u"m", z*u"m"], :free)
    n2 = Asap.Node([Lx*u"m", 0.0u"m", z*u"m"], :free)
    n3 = Asap.Node([Lx*u"m", Ly*u"m", z*u"m"], :free)
    n4 = Asap.Node([0.0u"m", Ly*u"m", z*u"m"], :free)
    # Center node
    nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", z*u"m"], :free)
    
    # 4 triangles meeting at center
    s1 = Asap.ShellTri3((n1, n2, nc), thickness*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
    s2 = Asap.ShellTri3((n2, n3, nc), thickness*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
    s3 = Asap.ShellTri3((n3, n4, nc), thickness*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
    s4 = Asap.ShellTri3((n4, n1, nc), thickness*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
    
    return [s1, s2, s3, s4], [n1, n2, n3, n4, nc]
end

"""Create beam elements around the perimeter of a rectangle."""
function make_rect_beams(Lx, Ly; z=0.0)
    n1 = Asap.Node([0.0u"m", 0.0u"m", z*u"m"], :free)
    n2 = Asap.Node([Lx*u"m", 0.0u"m", z*u"m"], :free)
    n3 = Asap.Node([Lx*u"m", Ly*u"m", z*u"m"], :free)
    n4 = Asap.Node([0.0u"m", Ly*u"m", z*u"m"], :free)
    
    sec = Asap.Section(
        0.005u"m^2", 200e9u"Pa", 80e9u"Pa",
        1.0e-4u"m^4", 5.0e-5u"m^4", 2.0e-7u"m^4", 7850u"kg/m^3"
    )
    
    beam1 = Asap.Element(n1, n2, sec, :beam_bottom)  # Edge 1: bottom
    beam2 = Asap.Element(n2, n3, sec, :beam_right)   # Edge 2: right  
    beam3 = Asap.Element(n3, n4, sec, :beam_top)     # Edge 3: top
    beam4 = Asap.Element(n4, n1, sec, :beam_left)    # Edge 4: left
    
    return [beam1, beam2, beam3, beam4], [n1, n2, n3, n4]
end

# =============================================================================
# 1. Edge Matching Tests
# =============================================================================

@testset "Edge Matching" begin
    
    @testset "Exact match - same orientation" begin
        # Polygon: square with vertices at (0,0), (4,0), (4,3), (0,3)
        vertices = [Point(0.0u"m", 0.0u"m"), Point(4.0u"m", 0.0u"m"),
                   Point(4.0u"m", 3.0u"m"), Point(0.0u"m", 3.0u"m")]
        
        # Beam from (0,0) to (4,0) should match edge 1
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 80e9u"Pa", 
                          1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4", 7850u"kg/m^3")
        beam = Asap.Element(n1, n2, sec, :test)
        
        edge_idx = Asap.match_beam_to_polygon_edge(beam, vertices)
        @test edge_idx == 1
    end
    
    @testset "Exact match - reversed orientation" begin
        vertices = [Point(0.0u"m", 0.0u"m"), Point(4.0u"m", 0.0u"m"),
                   Point(4.0u"m", 3.0u"m"), Point(0.0u"m", 3.0u"m")]
        
        # Beam from (4,0) to (0,0) - reversed direction
        n1 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4", 7850u"kg/m^3")
        beam = Asap.Element(n1, n2, sec, :test)
        
        edge_idx = Asap.match_beam_to_polygon_edge(beam, vertices)
        @test edge_idx == 1  # Should still match edge 1
    end
    
    @testset "Match all four edges" begin
        vertices = [Point(0.0u"m", 0.0u"m"), Point(4.0u"m", 0.0u"m"),
                   Point(4.0u"m", 3.0u"m"), Point(0.0u"m", 3.0u"m")]
        
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4", 7850u"kg/m^3")
        
        # Create beams for all 4 edges
        beams_data = [
            ((0.0, 0.0), (4.0, 0.0), 1),  # Bottom
            ((4.0, 0.0), (4.0, 3.0), 2),  # Right
            ((4.0, 3.0), (0.0, 3.0), 3),  # Top
            ((0.0, 3.0), (0.0, 0.0), 4),  # Left
        ]
        
        for (start, stop, expected_edge) in beams_data
            n1 = Asap.Node([start[1]*u"m", start[2]*u"m", 0.0u"m"], :free)
            n2 = Asap.Node([stop[1]*u"m", stop[2]*u"m", 0.0u"m"], :free)
            beam = Asap.Element(n1, n2, sec, :test)
            
            edge_idx = Asap.match_beam_to_polygon_edge(beam, vertices)
            @test edge_idx == expected_edge
        end
    end
    
    @testset "No match - interior beam" begin
        vertices = [Point(0.0u"m", 0.0u"m"), Point(4.0u"m", 0.0u"m"),
                   Point(4.0u"m", 3.0u"m"), Point(0.0u"m", 3.0u"m")]
        
        # Beam in the interior (not on any edge)
        n1 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([3.0u"m", 2.0u"m", 0.0u"m"], :free)
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4", 7850u"kg/m^3")
        beam = Asap.Element(n1, n2, sec, :interior)
        
        edge_idx = Asap.match_beam_to_polygon_edge(beam, vertices)
        @test isnothing(edge_idx)
    end
    
    @testset "Tolerance handling" begin
        vertices = [Point(0.0u"m", 0.0u"m"), Point(4.0u"m", 0.0u"m"),
                   Point(4.0u"m", 3.0u"m"), Point(0.0u"m", 3.0u"m")]
        
        sec = Asap.Section(0.01u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4", 7850u"kg/m^3")
        
        # Beam endpoints slightly off (1e-9 m)
        n1 = Asap.Node([1e-9*u"m", 1e-9*u"m", 0.0u"m"], :free)
        n2 = Asap.Node([(4.0 - 1e-9)*u"m", 1e-9*u"m", 0.0u"m"], :free)
        beam = Asap.Element(n1, n2, sec, :test)
        
        # Should match with default tolerance
        edge_idx = Asap.match_beam_to_polygon_edge(beam, vertices)
        @test edge_idx == 1
        
        # Should not match with very tight tolerance
        edge_idx_strict = Asap.match_beam_to_polygon_edge(beam, vertices; tol=1e-12)
        @test isnothing(edge_idx_strict)
    end
    
    @testset "Tuple-based API" begin
        vertices = [(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]
        
        # Test tuple-based edge matching
        edge_idx = Asap.match_beam_to_polygon_edge((0.0, 0.0), (4.0, 0.0), vertices)
        @test edge_idx == 1
        
        edge_idx = Asap.match_beam_to_polygon_edge((4.0, 0.0), (4.0, 3.0), vertices)
        @test edge_idx == 2
    end
    
end

# =============================================================================
# 2. Shell Boundary Extraction Tests  
# =============================================================================

@testset "Shell Boundary Extraction" begin
    
    @testset "Four triangles meeting at center" begin
        shells, nodes = make_rect_shells(4.0, 3.0)
        
        boundary = Asap.get_shell_panel_boundary(shells)
        
        @test length(boundary) == 4  # 4 vertices
        
        # Check area of boundary polygon (should be 4 × 3 = 12)
        area = Asap._polygon_area(boundary)
        @test isapprox(abs(area), 12.0, rtol=0.01)
    end
    
    @testset "Two-triangle mesh (diagonal split)" begin
        # Rectangle split diagonally: 2 triangles
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([4.0u"m", 3.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 3.0u"m", 0.0u"m"], :free)
        
        s1 = Asap.ShellTri3((n1, n2, n3), 0.15u"m", 30e9u"Pa", 0.2; id=:slab, ρ=2400.0)
        s2 = Asap.ShellTri3((n1, n3, n4), 0.15u"m", 30e9u"Pa", 0.2; id=:slab, ρ=2400.0)
        
        boundary = Asap.get_shell_panel_boundary([s1, s2])
        
        @test length(boundary) == 4
        area = Asap._polygon_area(boundary)
        @test isapprox(abs(area), 12.0, rtol=0.01)
    end
    
    @testset "Panel ID filtering" begin
        # Create shells with different IDs
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([2.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 2.0u"m", 0.0u"m"], :free)
        nc = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        # Panel A
        sa1 = Asap.ShellTri3((n1, n2, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_A, ρ=2400.0)
        sa2 = Asap.ShellTri3((n2, n3, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_A, ρ=2400.0)
        
        # Panel B (different id)
        sb1 = Asap.ShellTri3((n3, n4, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_B, ρ=2400.0)
        sb2 = Asap.ShellTri3((n4, n1, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_B, ρ=2400.0)
        
        all_shells = [sa1, sa2, sb1, sb2]
        
        # Filter by panel_A
        boundary_A = Asap.get_shell_panel_boundary(all_shells; panel_id=:panel_A)
        @test length(boundary_A) >= 3  # At least a triangle
        
        # Filter by panel_B
        boundary_B = Asap.get_shell_panel_boundary(all_shells; panel_id=:panel_B)
        @test length(boundary_B) >= 3
    end
    
end

# =============================================================================
# 3. Shell-to-Tributary Conversion Tests
# =============================================================================

@testset "Shell to Tributary Conversion" begin
    
    @testset "Simple rectangle - area conservation" begin
        Lx, Ly = 4.0, 3.0
        shells, shell_nodes = make_rect_shells(Lx, Ly)
        beams, beam_nodes = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure)
        
        @test length(loads) == 4  # One load per beam
        
        # Check total tributary area via load intensity integration
        # For uniform pressure, total force = pressure × area
        # Each TributaryLoad integrates to force on that beam
        total_force = 0.0
        for load in loads
            # Approximate integration using trapezoidal rule on intensities
            intensities_vals = Asap.intensities(load)
            # Compute beam length from node positions (element.length is 0 until processed)
            beam_length = beam_length_from_nodes(load.element)
            positions = load.positions
            
            # Trapezoidal integration
            for i in 1:(length(positions)-1)
                ds = (positions[i+1] - positions[i]) * beam_length
                avg_intensity = (intensities_vals[i] + intensities_vals[i+1]) / 2
                total_force += avg_intensity * ds
            end
        end
        
        expected_force = 5000.0 * Lx * Ly  # 60000 N
        @test isapprox(total_force, expected_force, rtol=0.02)
    end
    
    @testset "Tributary widths physically reasonable" begin
        Lx, Ly = 4.0, 3.0  # Short span = 3m
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure)
        
        # For isotropic straight skeleton of 4×3 rectangle:
        # Max tributary width = min(Lx, Ly) / 2 = 1.5 m
        expected_max_width = min(Lx, Ly) / 2
        
        for load in loads
            max_width = maximum(ustrip.(u"m", load.widths))
            @test max_width <= expected_max_width * 1.01  # Allow 1% tolerance
        end
    end
    
    @testset "TributaryLoad intensities" begin
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"  # 5 kPa
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure)
        
        # Max intensity = pressure × max_width = 5000 × 1.5 = 7500 N/m
        expected_max_intensity = 5000.0 * min(Lx, Ly) / 2
        
        for load in loads
            intensities_vals = Asap.intensities(load)
            max_intensity = maximum(intensities_vals)
            @test max_intensity <= expected_max_intensity * 1.01
        end
    end
    
    @testset "Unmatched beam throws error" begin
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        
        # Create beams that DON'T match the panel edges
        n1 = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)
        sec = Asap.Section(0.005u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 5e-5u"m^4", 2e-7u"m^4", 7850u"kg/m^3")
        interior_beam = Asap.Element(n1, n2, sec, :interior)
        
        pressure = 5000.0u"Pa"
        
        # Should throw an error
        @test_throws Exception Asap.shell_to_tributary_loads(shells, [interior_beam], pressure)
    end
    
end

# =============================================================================
# 4. One-Way vs Isotropic Distribution Tests
# =============================================================================

@testset "One-Way vs Isotropic Distribution" begin
    
    @testset "Two-way isotropic: all beams get load" begin
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure; axis=nothing)
        
        # All 4 beams should have non-zero load
        for load in loads
            max_intensity = maximum(Asap.intensities(load))
            @test max_intensity > 0.0
            @info "Isotropic: $(load.element.id) max intensity = $max_intensity N/m"
        end
    end
    
    @testset "One-way X: load only on Y-parallel edges" begin
        # Slab spans along X → load goes to edges perpendicular to X (Y-parallel edges)
        # Y-parallel edges: left (x=0) and right (x=Lx)
        # X-parallel edges: bottom (y=0) and top (y=Ly) → should get ZERO load
        
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure; axis=[1.0, 0.0])
        
        for load in loads
            max_intensity = maximum(Asap.intensities(load))
            beam_id = load.element.id
            
            if beam_id == :beam_bottom || beam_id == :beam_top
                # X-parallel edges: should get ZERO load
                @test max_intensity < 1.0  # Effectively zero (allow numerical noise)
                @info "One-way X: $beam_id (X-parallel) intensity = $max_intensity N/m ✓ (expected ~0)"
            else
                # Y-parallel edges: should get ALL load
                @test max_intensity > 1000.0  # Significant load
                @info "One-way X: $beam_id (Y-parallel) intensity = $max_intensity N/m ✓ (expected >0)"
            end
        end
    end
    
    @testset "One-way Y: load only on X-parallel edges" begin
        # Slab spans along Y → load goes to edges perpendicular to Y (X-parallel edges)
        # X-parallel edges: bottom (y=0) and top (y=Ly)
        # Y-parallel edges: left (x=0) and right (x=Lx) → should get ZERO load
        
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure; axis=[0.0, 1.0])
        
        for load in loads
            max_intensity = maximum(Asap.intensities(load))
            beam_id = load.element.id
            
            if beam_id == :beam_left || beam_id == :beam_right
                # Y-parallel edges: should get ZERO load
                @test max_intensity < 1.0  # Effectively zero
                @info "One-way Y: $beam_id (Y-parallel) intensity = $max_intensity N/m ✓ (expected ~0)"
            else
                # X-parallel edges: should get ALL load  
                @test max_intensity > 1000.0  # Significant load
                @info "One-way Y: $beam_id (X-parallel) intensity = $max_intensity N/m ✓ (expected >0)"
            end
        end
    end
    
    @testset "One-way vs isotropic: total load same" begin
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        pressure = 5000.0u"Pa"
        
        # Helper to compute total force from loads
        function total_force(loads)
            force = 0.0
            for load in loads
                intensities_vals = Asap.intensities(load)
                beam_length = beam_length_from_nodes(load.element)
                positions = load.positions
                for i in 1:(length(positions)-1)
                    ds = (positions[i+1] - positions[i]) * beam_length
                    avg = (intensities_vals[i] + intensities_vals[i+1]) / 2
                    force += avg * ds
                end
            end
            return force
        end
        
        loads_iso = Asap.shell_to_tributary_loads(shells, beams, pressure; axis=nothing)
        loads_x = Asap.shell_to_tributary_loads(shells, beams, pressure; axis=[1.0, 0.0])
        loads_y = Asap.shell_to_tributary_loads(shells, beams, pressure; axis=[0.0, 1.0])
        
        F_iso = total_force(loads_iso)
        F_x = total_force(loads_x)
        F_y = total_force(loads_y)
        
        expected = 5000.0 * Lx * Ly  # 60000 N
        
        # Isotropic distribution: total force = pressure × area (straight skeleton covers full area)
        @test isapprox(F_iso, expected, rtol=0.02)
        
        # One-way distribution: for rectangular panels, one-way X and Y both cover full area
        # - One-way X puts load on Y-parallel edges (each gets width = Lx/2)
        # - One-way Y puts load on X-parallel edges (each gets width = Ly/2)
        # Both should equal the total slab area when integrated
        @test isapprox(F_x, expected, rtol=0.02)  # Full load via Y-edges
        @test isapprox(F_y, expected, rtol=0.02)  # Full load via X-edges
        
        # For rectangular panels, F_x ≈ F_y ≈ expected
        @test isapprox(F_x, F_y, rtol=0.02)
        
        @info "Total force comparison" F_iso F_x F_y expected
    end
    
    @testset "Symbolic axis API" begin
        Lx, Ly = 4.0, 3.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        pressure = 5000.0u"Pa"
        
        # Test symbolic axis specification
        loads_iso = Asap.shell_to_tributary_loads(shells, beams, pressure, :isotropic)
        loads_x = Asap.shell_to_tributary_loads(shells, beams, pressure, :x)
        loads_y = Asap.shell_to_tributary_loads(shells, beams, pressure, :y)
        
        @test length(loads_iso) == 4
        @test length(loads_x) == 4
        @test length(loads_y) == 4
    end
    
end

# =============================================================================
# 5. Multi-Panel Support (Interior Beams)
# =============================================================================

@testset "Multi-Panel Support" begin
    
    @testset "Two panels sharing interior beam" begin
        # Two adjacent panels: Panel A (left) and Panel B (right)
        # They share an interior beam at x = 2
        #
        #   Panel A          Panel B
        #  (0,0)-(2,0)     (2,0)-(4,0)
        #    |     |         |     |
        #  (0,3)-(2,3)     (2,3)-(4,3)
        #
        # Interior beam: (2,0) to (2,3)
        
        # Panel A shells
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([2.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([2.0u"m", 3.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 3.0u"m", 0.0u"m"], :free)
        nc_a = Asap.Node([1.0u"m", 1.5u"m", 0.0u"m"], :free)
        
        sa1 = Asap.ShellTri3((n1, n2, nc_a), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_A, ρ=2400.0)
        sa2 = Asap.ShellTri3((n2, n3, nc_a), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_A, ρ=2400.0)
        sa3 = Asap.ShellTri3((n3, n4, nc_a), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_A, ρ=2400.0)
        sa4 = Asap.ShellTri3((n4, n1, nc_a), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_A, ρ=2400.0)
        panel_A = [sa1, sa2, sa3, sa4]
        
        # Panel B shells
        n5 = Asap.Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
        n6 = Asap.Node([4.0u"m", 3.0u"m", 0.0u"m"], :free)
        nc_b = Asap.Node([3.0u"m", 1.5u"m", 0.0u"m"], :free)
        
        # Reuse n2, n3 as shared nodes
        sb1 = Asap.ShellTri3((n2, n5, nc_b), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_B, ρ=2400.0)
        sb2 = Asap.ShellTri3((n5, n6, nc_b), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_B, ρ=2400.0)
        sb3 = Asap.ShellTri3((n6, n3, nc_b), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_B, ρ=2400.0)
        sb4 = Asap.ShellTri3((n3, n2, nc_b), 0.15u"m", 30e9u"Pa", 0.2; id=:panel_B, ρ=2400.0)
        panel_B = [sb1, sb2, sb3, sb4]
        
        # Create beams: perimeter + interior
        sec = Asap.Section(0.005u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 5e-5u"m^4", 2e-7u"m^4", 7850u"kg/m^3")
        
        # Exterior beams (panel A)
        beam_a_bottom = Asap.Element(n1, n2, sec, :a_bottom)
        beam_a_left = Asap.Element(n4, n1, sec, :a_left)
        beam_a_top = Asap.Element(n3, n4, sec, :a_top)
        
        # Interior beam (shared)
        beam_interior = Asap.Element(n2, n3, sec, :interior)
        
        # Exterior beams (panel B)
        beam_b_bottom = Asap.Element(n2, n5, sec, :b_bottom)
        beam_b_right = Asap.Element(n5, n6, sec, :b_right)
        beam_b_top = Asap.Element(n6, n3, sec, :b_top)
        
        all_beams = [beam_a_bottom, beam_a_left, beam_a_top, beam_interior,
                     beam_b_bottom, beam_b_right, beam_b_top]
        
        # Use multi-panel function
        pressure = 5000.0u"Pa"
        loads = Asap.shell_panels_to_tributary_loads([panel_A, panel_B], all_beams, pressure)
        
        @test length(loads) == 7  # All beams get a load
        
        # Find interior beam load
        interior_load = nothing
        for load in loads
            if load.element.id == :interior
                interior_load = load
                break
            end
        end
        
        @test !isnothing(interior_load)
        
        # Interior beam should have larger tributary width (sum from both sides)
        max_width_interior = maximum(ustrip.(u"m", interior_load.widths))
        
        # For 2m × 3m panels, max width from one side ≈ 1m
        # Interior beam gets from both sides, so max ≈ 2m
        @test max_width_interior > 1.5  # Should be significantly more than single-panel
        @info "Interior beam max tributary width: $max_width_interior m (expected ~2m)"
        
        # Total load should equal total area × pressure
        total_area = 2 * 2.0 * 3.0  # Two 2×3 panels = 12 m²
        expected_force = 5000.0 * total_area
        
        actual_force = 0.0
        for load in loads
            intensities_vals = Asap.intensities(load)
            beam_length = beam_length_from_nodes(load.element)
            positions = load.positions
            for i in 1:(length(positions)-1)
                ds = (positions[i+1] - positions[i]) * beam_length
                avg = (intensities_vals[i] + intensities_vals[i+1]) / 2
                actual_force += avg * ds
            end
        end
        
        @test isapprox(actual_force, expected_force, rtol=0.05)
        @info "Multi-panel total force: $actual_force N (expected: $expected_force N)"
    end
    
    @testset "Unmatched beam in multi-panel throws error" begin
        # Create two panels but include a beam that doesn't match either
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([2.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)
        n4 = Asap.Node([0.0u"m", 2.0u"m", 0.0u"m"], :free)
        nc = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        s1 = Asap.ShellTri3((n1, n2, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel, ρ=2400.0)
        s2 = Asap.ShellTri3((n2, n3, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel, ρ=2400.0)
        s3 = Asap.ShellTri3((n3, n4, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel, ρ=2400.0)
        s4 = Asap.ShellTri3((n4, n1, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:panel, ρ=2400.0)
        panel = [s1, s2, s3, s4]
        
        # Unmatched beam (in the middle, not on boundary)
        n_mid1 = Asap.Node([0.5u"m", 0.5u"m", 0.0u"m"], :free)
        n_mid2 = Asap.Node([1.5u"m", 1.5u"m", 0.0u"m"], :free)
        sec = Asap.Section(0.005u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 5e-5u"m^4", 2e-7u"m^4", 7850u"kg/m^3")
        unmatched_beam = Asap.Element(n_mid1, n_mid2, sec, :unmatched)
        
        pressure = 5000.0u"Pa"
        
        @test_throws Exception Asap.shell_panels_to_tributary_loads([panel], [unmatched_beam], pressure)
    end
    
end

# =============================================================================
# 6. Global Equilibrium Tests
# =============================================================================

@testset "Global Equilibrium" begin
    
    @testset "Tributary model equilibrium" begin
        # Create a complete structural model with columns, beams, and tributary loads
        # Verify total reactions = total applied load
        
        Lx, Ly, H = 4.0, 3.0, 3.0
        pressure = 5000.0  # Pa
        
        # Nodes - column bases (fixed)
        n1_base = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2_base = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3_base = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4_base = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        # Nodes - beam corners (free)
        n1_top = Asap.Node([0.0u"m", 0.0u"m", H*u"m"], :free)
        n2_top = Asap.Node([Lx*u"m", 0.0u"m", H*u"m"], :free)
        n3_top = Asap.Node([Lx*u"m", Ly*u"m", H*u"m"], :free)
        n4_top = Asap.Node([0.0u"m", Ly*u"m", H*u"m"], :free)
        
        # Sections
        col_sec = Asap.Section(6.25e-3u"m^2", 200e9u"Pa", 80e9u"Pa",
                              2.0e-4u"m^4", 7.0e-5u"m^4", 3.0e-7u"m^4", 7850u"kg/m^3")
        beam_sec = Asap.Section(4.94e-3u"m^2", 200e9u"Pa", 80e9u"Pa",
                               1.5e-4u"m^4", 5.0e-5u"m^4", 2.0e-7u"m^4", 7850u"kg/m^3")
        
        # Columns
        col1 = Asap.Element(n1_base, n1_top, col_sec, :col1)
        col2 = Asap.Element(n2_base, n2_top, col_sec, :col2)
        col3 = Asap.Element(n3_base, n3_top, col_sec, :col3)
        col4 = Asap.Element(n4_base, n4_top, col_sec, :col4)
        
        # Beams
        beam1 = Asap.Element(n1_top, n2_top, beam_sec, :beam_bottom)
        beam2 = Asap.Element(n2_top, n3_top, beam_sec, :beam_right)
        beam3 = Asap.Element(n3_top, n4_top, beam_sec, :beam_top)
        beam4 = Asap.Element(n4_top, n1_top, beam_sec, :beam_left)
        
        # Create shells for tributary computation
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", H*u"m"], :free)
        shells = [
            Asap.ShellTri3((n1_top, n2_top, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:slab, ρ=2400.0),
            Asap.ShellTri3((n2_top, n3_top, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:slab, ρ=2400.0),
            Asap.ShellTri3((n3_top, n4_top, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:slab, ρ=2400.0),
            Asap.ShellTri3((n4_top, n1_top, nc), 0.15u"m", 30e9u"Pa", 0.2; id=:slab, ρ=2400.0),
        ]
        beams = [beam1, beam2, beam3, beam4]
        
        # Generate tributary loads
        trib_loads = Asap.shell_to_tributary_loads(shells, beams, pressure*u"Pa")
        
        # Build and solve model
        nodes = [n1_base, n2_base, n3_base, n4_base, n1_top, n2_top, n3_top, n4_top]
        elements = [col1, col2, col3, col4, beam1, beam2, beam3, beam4]
        
        model = Asap.FrameModel(nodes, elements, trib_loads)
        Asap.process!(model)
        Asap.solve!(model)
        Asap.post_process!(model)
        Asap.reactions!(model)
        
        # Check equilibrium
        total_load = pressure * Lx * Ly
        R_z = sum(ustrip(u"N", n.reaction[3]) for n in [n1_base, n2_base, n3_base, n4_base])
        
        @test isapprox(R_z, total_load, rtol=0.01)
        @info "Equilibrium check: R_z = $R_z N, expected = $total_load N"
        
        # Horizontal reactions should be near zero (vertical load only)
        R_x = sum(ustrip(u"N", n.reaction[1]) for n in [n1_base, n2_base, n3_base, n4_base])
        R_y = sum(ustrip(u"N", n.reaction[2]) for n in [n1_base, n2_base, n3_base, n4_base])
        
        @test abs(R_x) < total_load * 0.001  # < 0.1% of vertical load
        @test abs(R_y) < total_load * 0.001
    end
    
end

# =============================================================================
# 7. Edge Cases
# =============================================================================

@testset "Edge Cases" begin
    
    @testset "Empty shell list" begin
        beams, _ = make_rect_beams(4.0, 3.0)
        pressure = 5000.0u"Pa"
        
        loads = Asap.shell_to_tributary_loads(Asap.ShellTri3[], beams, pressure)
        @test isempty(loads)
    end
    
    @testset "Empty beam list" begin
        shells, _ = make_rect_shells(4.0, 3.0)
        pressure = 5000.0u"Pa"
        
        loads = Asap.shell_to_tributary_loads(shells, Asap.Element[], pressure)
        @test isempty(loads)
    end
    
    @testset "Single triangle panel" begin
        # Minimal case: 1 ShellTri3, 3 edge beams
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([3.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.5u"m", 2.0u"m", 0.0u"m"], :free)
        
        shell = Asap.ShellTri3((n1, n2, n3), 0.15u"m", 30e9u"Pa", 0.2; id=:tri, ρ=2400.0)
        
        sec = Asap.Section(0.005u"m^2", 200e9u"Pa", 80e9u"Pa",
                          1e-4u"m^4", 5e-5u"m^4", 2e-7u"m^4", 7850u"kg/m^3")
        
        beam1 = Asap.Element(n1, n2, sec, :b1)
        beam2 = Asap.Element(n2, n3, sec, :b2)
        beam3 = Asap.Element(n3, n1, sec, :b3)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads([shell], [beam1, beam2, beam3], pressure)
        
        @test length(loads) == 3
        
        # Total force should equal triangle area × pressure
        # Area = 0.5 × 3 × 2 = 3 m²
        triangle_area = 0.5 * 3.0 * 2.0
        expected_force = 5000.0 * triangle_area
        
        actual_force = 0.0
        for load in loads
            intensities_vals = Asap.intensities(load)
            beam_length = beam_length_from_nodes(load.element)
            positions = load.positions
            for i in 1:(length(positions)-1)
                ds = (positions[i+1] - positions[i]) * beam_length
                avg = (intensities_vals[i] + intensities_vals[i+1]) / 2
                actual_force += avg * ds
            end
        end
        
        @test isapprox(actual_force, expected_force, rtol=0.05)
    end
    
    @testset "High aspect ratio panel" begin
        # Very thin panel: 10m × 1m (aspect ratio 10:1)
        Lx, Ly = 10.0, 1.0
        shells, _ = make_rect_shells(Lx, Ly)
        beams, _ = make_rect_beams(Lx, Ly)
        
        pressure = 5000.0u"Pa"
        loads = Asap.shell_to_tributary_loads(shells, beams, pressure)
        
        @test length(loads) == 4
        
        # Max tributary width should be min(Lx, Ly)/2 = 0.5m
        for load in loads
            max_width = maximum(ustrip.(u"m", load.widths))
            @test max_width <= min(Lx, Ly) / 2 * 1.01
        end
    end
    
end
