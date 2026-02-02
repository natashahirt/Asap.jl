#=
Test: Shell Spatial Queries and Region Integration
==================================================
Tests the shell_tris_at_point, shell_tris_in_region, and
bending_moments dispatch extensions for polygon/point queries.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "Shell Queries & Region Integration" begin
    
    # Common setup: 4x3m plate with 4 triangular shells
    E = 30e9  # Pa (concrete)
    ν = 0.2
    ρ = 2400.0  # kg/m³
    t = 0.2  # m thickness
    
    Lx, Ly = 4.0, 3.0  # m
    pressure = 5000.0  # Pa (5 kPa)
    
    # Helper to create a standard test plate with denser mesh
    function create_test_plate(n_subdivisions=4)
        section = Asap.ShellSection(t*u"m", E*u"Pa", ν; ρ=ρ*u"kg/m^3")
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :pinned)
        n4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :pinned)
        
        shells = Asap.Shell((n1, n2, n3, n4), section; n=n_subdivisions)
        nodes = Asap.get_nodes(shells)
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        model = Asap.ShellModel(nodes, shells, [load])
        Asap.process!(model)
        Asap.solve!(model)
        
        return model, shells, nodes
    end
    
    @testset "shell_centroid" begin
        # Simple triangle
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([3.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([1.5u"m", 3.0u"m", 0.0u"m"], :free)
        
        shell = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        
        cx, cy = Asap.shell_centroid(shell)
        
        # Centroid of triangle with vertices (0,0), (3,0), (1.5,3) is (1.5, 1.0)
        @test isapprox(cx, 1.5, atol=1e-10)
        @test isapprox(cy, 1.0, atol=1e-10)
        
        @info "shell_centroid test passed" centroid=(cx, cy)
    end
    
    @testset "shell_centroid_3d" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 5.0u"m"], :free)
        n2 = Asap.Node([3.0u"m", 0.0u"m", 5.0u"m"], :free)
        n3 = Asap.Node([1.5u"m", 3.0u"m", 5.0u"m"], :free)
        
        shell = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        
        cx, cy, cz = Asap.shell_centroid_3d(shell)
        
        @test isapprox(cx, 1.5, atol=1e-10)
        @test isapprox(cy, 1.0, atol=1e-10)
        @test isapprox(cz, 5.0, atol=1e-10)
        
        @info "shell_centroid_3d test passed" centroid=(cx, cy, cz)
    end
    
    @testset "shell_tris_at_point - interior" begin
        model, shells, _ = create_test_plate(4)
        
        # Query at center of plate - should find at least 1 triangle
        center_pt = (Lx/2, Ly/2)
        tris_at_center = Asap.shell_tris_at_point(model, center_pt; tol=0.01)
        
        @test !isempty(tris_at_center)
        @test all(t isa Asap.ShellTri3 for t in tris_at_center)
        
        @info "shell_tris_at_point (interior) test passed" n_tris=length(tris_at_center)
    end
    
    @testset "shell_tris_at_point - edge (multiple triangles)" begin
        model, shells, _ = create_test_plate(4)
        
        # Query at a point that's likely on an interior edge between triangles
        # With n=4 subdivision, there should be internal edges
        edge_pt = (Lx/4, Ly/4)  # Quarter point
        tris_at_edge = Asap.shell_tris_at_point(model, edge_pt; tol=0.05)
        
        # At an edge, we might get 2 triangles (or 1 if not exactly on edge)
        @test !isempty(tris_at_edge)
        
        @info "shell_tris_at_point (edge) test passed" n_tris=length(tris_at_edge)
    end
    
    @testset "shell_tris_at_point - outside returns empty" begin
        model, shells, _ = create_test_plate(4)
        
        # Query outside the plate
        outside_pt = (Lx + 1.0, Ly + 1.0)
        tris_outside = Asap.shell_tris_at_point(model, outside_pt)
        
        @test isempty(tris_outside)
        
        @info "shell_tris_at_point (outside) test passed"
    end
    
    @testset "shell_tris_in_region - full plate" begin
        model, shells, _ = create_test_plate(4)
        
        # Region covering the entire plate
        full_region = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
        tris_in_full = Asap.shell_tris_in_region(model, full_region)
        
        # Should find all triangles
        @test length(tris_in_full) == length(shells)
        
        @info "shell_tris_in_region (full) test passed" n_tris=length(tris_in_full)
    end
    
    @testset "shell_tris_in_region - partial (strip)" begin
        model, shells, _ = create_test_plate(4)
        
        # Horizontal strip covering left half
        strip_region = [(0.0, 0.0), (Lx/2, 0.0), (Lx/2, Ly), (0.0, Ly)]
        tris_in_strip = Asap.shell_tris_in_region(model, strip_region)
        
        # Should find some but not all
        @test !isempty(tris_in_strip)
        @test length(tris_in_strip) < length(shells)
        
        # All found triangles should have centroid in the strip
        for tri in tris_in_strip
            cx, cy = Asap.shell_centroid(tri)
            @test 0.0 <= cx <= Lx/2 + 0.01  # Small tolerance
            @test 0.0 <= cy <= Ly + 0.01
        end
        
        @info "shell_tris_in_region (strip) test passed" n_tris=length(tris_in_strip)
    end
    
    @testset "bending_moments(elem, model)" begin
        model, shells, _ = create_test_plate(4)
        
        # Get moments for first shell element
        M = Asap.bending_moments(shells[1], model)
        
        @test length(M) == 3  # [Mxx, Myy, Mxy]
        @test all(isfinite.(M))
        
        @info "bending_moments(elem, model) test passed" M=M
    end
    
    @testset "bending_moments(Vector{ShellElement}, model)" begin
        model, shells, _ = create_test_plate(4)
        
        # Get moments for subset of shells
        subset = shells[1:min(5, length(shells))]
        Ms = Asap.bending_moments(subset, model)
        
        @test length(Ms) == length(subset)
        @test all(length(m) == 3 for m in Ms)
        
        @info "bending_moments(Vector, model) test passed" n_results=length(Ms)
    end
    
    @testset "bending_moments(model, pt) - single point" begin
        model, shells, _ = create_test_plate(4)
        
        # Query at plate center
        center_pt = (Lx/2, Ly/2)
        M = Asap.bending_moments(model, center_pt; tol=0.1)
        
        @test M !== nothing
        @test length(M) == 3
        @test all(isfinite.(M))
        
        @info "bending_moments(model, pt) test passed" M=M
    end
    
    @testset "bending_moments(model, pt) - outside returns nothing" begin
        model, shells, _ = create_test_plate(4)
        
        outside_pt = (Lx + 10.0, Ly + 10.0)
        M = Asap.bending_moments(model, outside_pt)
        
        @test M === nothing
        
        @info "bending_moments(model, pt) outside test passed"
    end
    
    @testset "bending_moments(model; polygon=...)" begin
        model, shells, _ = create_test_plate(4)
        
        # Full plate polygon
        full_poly = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
        result = Asap.bending_moments(model; polygon=full_poly)
        
        @test haskey(result, :Mxx)
        @test haskey(result, :Myy)
        @test haskey(result, :Mxy)
        @test haskey(result, :Mxx_avg)
        @test haskey(result, :Mxx_max)
        @test haskey(result, :area)
        @test haskey(result, :shell_tris)
        
        @test result.area > 0
        @test length(result.shell_tris) == length(shells)
        
        @info "bending_moments(model; polygon) test passed" area=result.area Mxx_max=result.Mxx_max
    end
    
    @testset "bending_moments(model; pts=...)" begin
        model, shells, _ = create_test_plate(4)
        
        # Multiple query points
        pts = [
            (Lx/4, Ly/4),
            (Lx/2, Ly/2),
            (3*Lx/4, 3*Ly/4),
            (Lx + 10.0, Ly + 10.0)  # Outside - should return nothing
        ]
        
        results = Asap.bending_moments(model; pts=pts, tol=0.1)
        
        @test length(results) == length(pts)
        # First 3 should have values
        @test results[1] !== nothing
        @test results[2] !== nothing
        @test results[3] !== nothing
        # Last one (outside) should be nothing
        @test results[4] === nothing
        
        @info "bending_moments(model; pts) test passed" n_results=count(r -> r !== nothing, results)
    end
    
    @testset "bending_moments(elems, model; polygon=...)" begin
        model, shells, _ = create_test_plate(4)
        
        # First get a subset of shells (left half)
        left_poly = [(0.0, 0.0), (Lx/2, 0.0), (Lx/2, Ly), (0.0, Ly)]
        left_tris = Asap.shell_tris_in_region(model, left_poly)
        
        # Now integrate over an even smaller region within that subset
        small_poly = [(0.0, 0.0), (Lx/4, 0.0), (Lx/4, Ly/2), (0.0, Ly/2)]
        result = Asap.bending_moments(left_tris, model; polygon=small_poly)
        
        @test result.area >= 0
        # The number of shells should be <= number in left_tris
        @test length(result.shell_tris) <= length(left_tris)
        
        @info "bending_moments(elems, model; polygon) test passed" n_in_region=length(result.shell_tris)
    end
    
    @testset "Region integration - conservation check" begin
        model, shells, _ = create_test_plate(4)
        
        # Split plate into left and right halves
        left_poly = [(0.0, 0.0), (Lx/2, 0.0), (Lx/2, Ly), (0.0, Ly)]
        right_poly = [(Lx/2, 0.0), (Lx, 0.0), (Lx, Ly), (Lx/2, Ly)]
        full_poly = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
        
        result_left = Asap.bending_moments(model; polygon=left_poly)
        result_right = Asap.bending_moments(model; polygon=right_poly)
        result_full = Asap.bending_moments(model; polygon=full_poly)
        
        # Total area should approximately equal sum of halves
        # (might not be exact due to triangles straddling the boundary)
        area_sum = result_left.area + result_right.area
        @test isapprox(area_sum, result_full.area, rtol=0.1)
        
        @info "Region integration conservation check passed" area_sum=area_sum area_full=result_full.area
    end
    
    @testset "Loaded plate - moments are non-trivial" begin
        model, shells, _ = create_test_plate(6)  # Finer mesh
        
        # Under uniform pressure, the plate should have non-zero bending moments
        center_pt = (Lx/2, Ly/2)
        M = Asap.bending_moments(model, center_pt; tol=0.1)
        
        @test M !== nothing
        # Under pressure loading, there should be significant bending
        # (values depend on formulation but should be non-zero)
        max_M = maximum(abs.(M))
        @test max_M > 0
        
        # Full plate integration
        full_poly = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
        result = Asap.bending_moments(model; polygon=full_poly)
        
        # Peak moments should be positive (in absolute value)
        @test result.Mxx_max > 0
        @test result.Myy_max > 0
        
        @info "Loaded plate moments test passed" M_center=M Mxx_max=result.Mxx_max
    end
end

println("\n✓ All Shell Queries & Region Integration tests passed!")
