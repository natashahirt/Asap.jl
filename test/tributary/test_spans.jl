# Tests for polygon span calculations

using Test
using Meshes
using Unitful
using Asap: get_polygon_span, SpanInfo, governing_spans, short_span, long_span, two_way_span

@testset "Polygon Span Calculations" begin
    
    # Helper to create Meshes.Point from (x, y) in meters
    point2d(x, y) = Point(x * u"m", y * u"m")
    
    @testset "Rectangle (axis-aligned)" begin
        # 10m x 6m rectangle
        verts = [point2d(0, 0), point2d(10, 0), point2d(10, 6), point2d(0, 6)]
        
        # One-way spans should match rectangle dimensions
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        
        @test span_x ≈ 10.0 atol=1e-9
        @test span_y ≈ 6.0 atol=1e-9
        
        # Two-way span = diagonal
        span_2way = get_polygon_span(verts)
        expected_diag = sqrt(10^2 + 6^2)
        @test span_2way ≈ expected_diag atol=1e-9
    end
    
    @testset "Square" begin
        # 8m x 8m square
        verts = [point2d(0, 0), point2d(8, 0), point2d(8, 8), point2d(0, 8)]
        
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        
        @test span_x ≈ 8.0 atol=1e-9
        @test span_y ≈ 8.0 atol=1e-9
        
        # Two-way span = diagonal = 8√2
        span_2way = get_polygon_span(verts)
        @test span_2way ≈ 8 * sqrt(2) atol=1e-9
    end
    
    @testset "Rotated Rectangle (45°)" begin
        # Rectangle rotated 45°, vertices at diamond orientation
        # Side length = 4√2, so width along original axes = 8
        s = 4 * sqrt(2)
        verts = [point2d(4, 0), point2d(8, 4), point2d(4, 8), point2d(0, 4)]
        
        # One-way spans along global axes
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        
        # Both should be 8 (bounding box width)
        @test span_x ≈ 8.0 atol=1e-6
        @test span_y ≈ 8.0 atol=1e-6
        
        # One-way span along diagonal (45°) should be the rectangle's actual width
        span_diag = get_polygon_span(verts; axis=[1.0, 1.0])
        @test span_diag ≈ s atol=1e-6
    end
    
    @testset "Parallelogram (skewed)" begin
        # Parallelogram: base 10m, height 5m, skewed 2m
        # (0,0) -> (10,0) -> (12,5) -> (2,5)
        verts = [point2d(0, 0), point2d(10, 0), point2d(12, 5), point2d(2, 5)]
        
        # Y-span should be the height = 5m (perpendicular distance)
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        @test span_y ≈ 5.0 atol=1e-6
        
        # X-span: from any vertex, max horizontal chord
        # From (0,0) going right hits (10,0) edge, then continues...
        # Actually the x-extent from leftmost to rightmost is 12m
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        @test span_x >= 10.0  # At least the base length
        @test span_x <= 12.0  # At most bounding box width
    end
    
    @testset "Trapezoid" begin
        # Isoceles trapezoid: bottom 10m, top 6m, height 4m
        # (0,0) -> (10,0) -> (8,4) -> (2,4)
        verts = [point2d(0, 0), point2d(10, 0), point2d(8, 4), point2d(2, 4)]
        
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        @test span_y ≈ 4.0 atol=1e-6
        
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        @test span_x ≈ 10.0 atol=1e-6  # Bottom edge is longest horizontal
        
        # Two-way: max vertex distance = diagonal from (0,0) to (8,4) or (10,0) to (2,4)
        span_2way = get_polygon_span(verts)
        diag1 = sqrt(8^2 + 4^2)  # (0,0) to (8,4)
        diag2 = sqrt(10^2 + 0^2)  # (0,0) to (10,0)
        diag3 = sqrt((10-2)^2 + 4^2)  # (10,0) to (2,4)
        @test span_2way ≈ max(diag1, diag2, diag3) atol=1e-6
    end
    
    @testset "L-Shape (non-convex)" begin
        # L-shape: 
        # (0,0) -> (6,0) -> (6,3) -> (3,3) -> (3,6) -> (0,6)
        verts = [point2d(0, 0), point2d(6, 0), point2d(6, 3), 
                 point2d(3, 3), point2d(3, 6), point2d(0, 6)]
        
        # Y-span along vertical: 6m (from bottom to top)
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        @test span_y ≈ 6.0 atol=1e-6
        
        # X-span along horizontal: 6m (from left to right)
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        @test span_x ≈ 6.0 atol=1e-6
        
        # Two-way: max vertex distance = diagonal from (0,0) to (6,3) or (0,6) to (6,0)
        span_2way = get_polygon_span(verts)
        diag1 = sqrt(6^2 + 3^2)  # (0,0) to (6,3)
        diag2 = sqrt(6^2 + 6^2)  # (0,0) to (6,6) - but (6,6) not a vertex
        diag3 = sqrt(6^2 + 6^2)  # (0,6) to (6,0)
        @test span_2way ≈ diag3 atol=1e-6  # sqrt(72)
    end
    
    @testset "Triangle" begin
        # Right triangle: (0,0), (8,0), (0,6)
        verts = [point2d(0, 0), point2d(8, 0), point2d(0, 6)]
        
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        
        @test span_x ≈ 8.0 atol=1e-6  # Base
        @test span_y ≈ 6.0 atol=1e-6  # Height
        
        # Two-way: hypotenuse = sqrt(8² + 6²) = 10
        span_2way = get_polygon_span(verts)
        @test span_2way ≈ 10.0 atol=1e-6
    end
    
    @testset "Offset Rectangle (translation invariance)" begin
        # Same 10x6 rectangle, but offset to (100, 200)
        verts = [point2d(100, 200), point2d(110, 200), 
                 point2d(110, 206), point2d(100, 206)]
        
        span_x = get_polygon_span(verts; axis=[1.0, 0.0])
        span_y = get_polygon_span(verts; axis=[0.0, 1.0])
        
        @test span_x ≈ 10.0 atol=1e-9
        @test span_y ≈ 6.0 atol=1e-9
    end
    
    @testset "Edge cases" begin
        # Single point
        verts_1 = [point2d(0, 0)]
        @test get_polygon_span(verts_1) == 0.0
        @test get_polygon_span(verts_1; axis=[1.0, 0.0]) == 0.0
        
        # Two points (line segment)
        verts_2 = [point2d(0, 0), point2d(5, 0)]
        @test get_polygon_span(verts_2) ≈ 5.0 atol=1e-9
        @test get_polygon_span(verts_2; axis=[1.0, 0.0]) == 0.0  # Not enough vertices for polygon
        
        # Zero-length axis
        verts_rect = [point2d(0, 0), point2d(10, 0), point2d(10, 6), point2d(0, 6)]
        span_zero = get_polygon_span(verts_rect; axis=[0.0, 0.0])
        # Should fall back to two-way
        @test span_zero ≈ sqrt(10^2 + 6^2) atol=1e-9
    end
    
    # =========================================================================
    # SpanInfo Tests
    # =========================================================================
    
    @testset "SpanInfo - auto axis detection" begin
        # 10m x 6m rectangle → primary should be 6m (short), secondary 10m (long)
        verts = [point2d(0, 0), point2d(10, 0), point2d(10, 6), point2d(0, 6)]
        si = SpanInfo(verts)
        
        @test si.primary ≈ 6.0 atol=1e-9
        @test si.secondary ≈ 10.0 atol=1e-9
        @test si.axis == (0.0, 1.0)  # Y-axis is short direction
        @test si.isotropic ≈ sqrt(10^2 + 6^2) atol=1e-9
        
        # Accessors
        @test short_span(si) ≈ 6.0 atol=1e-9
        @test long_span(si) ≈ 10.0 atol=1e-9
        @test two_way_span(si) ≈ sqrt(136) atol=1e-9
    end
    
    @testset "SpanInfo - custom axis" begin
        # 10m x 6m rectangle with custom axis along X
        verts = [point2d(0, 0), point2d(10, 0), point2d(10, 6), point2d(0, 6)]
        si = SpanInfo(verts; axis=(1.0, 0.0))
        
        @test si.primary ≈ 10.0 atol=1e-9   # X span
        @test si.secondary ≈ 6.0 atol=1e-9   # Y span (perpendicular)
        @test si.axis == (1.0, 0.0)
    end
    
    @testset "SpanInfo - 45° axis" begin
        # Square with 45° axis
        verts = [point2d(0, 0), point2d(8, 0), point2d(8, 8), point2d(0, 8)]
        si = SpanInfo(verts; axis=(1.0, 1.0))  # 45° diagonal
        
        # Span along diagonal should be 8√2
        @test si.primary ≈ 8 * sqrt(2) atol=1e-6
        @test si.secondary ≈ 8 * sqrt(2) atol=1e-6  # perpendicular diagonal is same for square
    end
    
    @testset "SpanInfo - governing_spans" begin
        # Three cells with different spans
        verts1 = [point2d(0, 0), point2d(8, 0), point2d(8, 6), point2d(0, 6)]
        verts2 = [point2d(0, 0), point2d(10, 0), point2d(10, 5), point2d(0, 5)]
        verts3 = [point2d(0, 0), point2d(7, 0), point2d(7, 9), point2d(0, 9)]
        
        si1 = SpanInfo(verts1)  # 6m x 8m → primary=6, secondary=8
        si2 = SpanInfo(verts2)  # 5m x 10m → primary=5, secondary=10
        si3 = SpanInfo(verts3)  # 7m x 9m → primary=7, secondary=9
        
        gov = governing_spans([si1, si2, si3])
        
        # Governing = max of each component
        @test gov.primary ≈ 7.0 atol=1e-9   # max(6, 5, 7)
        @test gov.secondary ≈ 10.0 atol=1e-9  # max(8, 10, 9)
        @test gov.isotropic ≈ sqrt(7^2 + 9^2) atol=1e-6  # max diagonal
        @test gov.axis == si1.axis  # Uses first cell's axis
    end
    
    @testset "SpanInfo - single cell governing" begin
        verts = [point2d(0, 0), point2d(5, 0), point2d(5, 5), point2d(0, 5)]
        si = SpanInfo(verts)
        
        gov = governing_spans([si])
        @test gov.primary == si.primary
        @test gov.secondary == si.secondary
        @test gov.isotropic == si.isotropic
    end
    
    println("✓ All span calculation tests passed!")
end
