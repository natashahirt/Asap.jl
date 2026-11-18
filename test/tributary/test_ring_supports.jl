# A beam ring is a supported hole. An opening of the same outline is not.
# One-way load splits at the midline between a ring edge and the outer beam.
# The opening takes none of that strip; the outer beam takes all of it.

using Test
using Asap
using Meshes: Point
using Unitful: @u_str

_pts(xy) = [Point(x * u"m", y * u"m") for (x, y) in xy]

function _partition(hole; support::Bool)
    outer = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    return get_tributary_polygons(_pts(outer); axis = [1.0, 0.0], holes = [_pts(hole)],
                                  hole_supports = [support])
end

_area(tribs, pred) = sum(t.area for t in tribs if pred(t.local_edge_idx); init = 0.0)

@testset "Beam ring supports" begin
    small = [(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0)]
    big = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)]

    ring = _partition(small; support = true)
    @test _area(ring, _ -> true) ≈ 96.0 atol = 1e-6
    # Vertical ring edges, 2 m tall, each takes half of the 4 m gap: 4 m² each.
    @test _area(ring, i -> i > 4) ≈ 8.0 atol = 1e-6
    again = _partition(small; support = true)
    @test _area(again, i -> i > 4) ≈ _area(ring, i -> i > 4)

    # A taller ring intercepts more of the bay, so it carries more slab.
    larger = _partition(big; support = true)
    @test _area(larger, i -> i > 4) ≈ 12.0 atol = 1e-6
    @test _area(larger, i -> i > 4) > _area(ring, i -> i > 4) + 1.0

    opening = _partition(small; support = false)
    @test _area(opening, i -> i > 4) ≈ 0.0 atol = 1e-9
    @test _area(opening, _ -> true) ≈ 96.0 atol = 1e-6
    # The same void, with no ring to share it, stays on the outer beams.
    @test _area(opening, i -> i <= 4) > _area(ring, i -> i <= 4) + 1.0
end

function _iso(hole; support::Bool)
    outer = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    return get_tributary_polygons(_pts(outer); holes = [_pts(hole)],
                                  hole_supports = [support])
end

@testset "Two-way beam ring supports" begin
    # Straight skeleton of a square frame. Gap g, outer edge
    # (2L − g) g / 4, hole edge the trapezoid between the ring and the midline.
    small = [(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0)]
    big = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)]

    ring = _iso(small; support = true)
    areas = Dict(t.local_edge_idx => t.area for t in ring)
    @test areas[1] ≈ 16.0 atol = 1e-4
    @test areas[5] ≈ 8.0 atol = 1e-4
    @test _area(ring, i -> i > 4) ≈ 32.0 atol = 1e-4
    @test _area(ring, _ -> true) ≈ 96.0 atol = 1e-4
    again = _iso(small; support = true)
    @test _area(again, i -> i > 4) ≈ _area(ring, i -> i > 4)

    larger = _iso(big; support = true)
    @test _area(larger, i -> i > 4) ≈ 33.0 atol = 1e-4
    @test _area(larger, i -> i > 4) > _area(ring, i -> i > 4)

    # The ring sits against the left side, so that beam takes less slab than
    # the far one.
    shifted = _iso([(1.0, 4.0), (3.0, 4.0), (3.0, 6.0), (1.0, 6.0)]; support = true)
    shift_a = Dict(t.local_edge_idx => t.area for t in shifted)
    @test shift_a[4] < shift_a[2]
    @test _area(shifted, i -> i > 4) > 0
    @test _area(shifted, _ -> true) ≈ 96.0 atol = 1e-4

    # A diamond is not parallel to the slab edge. The ring still takes slab,
    # and a second call matches. Its area is not the square-frame trapezoid.
    diamond = [(5.0, 3.0), (7.0, 5.0), (5.0, 7.0), (3.0, 5.0)]
    turned = _iso(diamond; support = true)
    @test _area(turned, _ -> true) ≈ 92.0 atol = 1e-4
    @test _area(turned, i -> i > 4) > 1.0
    @test _area(_iso(diamond; support = true), i -> i > 4) ≈ _area(turned, i -> i > 4)

    opening = _iso(small; support = false)
    @test _area(opening, i -> i > 4) ≈ 0.0 atol = 1e-9
    @test _area(opening, _ -> true) ≈ 96.0 atol = 1e-6
    @test _area(opening, i -> i <= 4) > _area(ring, i -> i <= 4)

    # The slab edge is split, and the ring edge overhangs that joint on
    # both pieces. Neither edge contains the other. The slab outside the
    # hole is still the whole tributary: 10×10 minus the 6×2 hole.
    split_outer = [(0.0, 0.0), (6.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    across = [(2.0, 4.0), (8.0, 4.0), (8.0, 6.0), (2.0, 6.0)]
    crossed = get_tributary_polygons(_pts(split_outer); holes = [_pts(across)],
                                     hole_supports = [true])
    @test _area(crossed, _ -> true) ≈ 88.0 atol = 1e-3
    @test _area(crossed, i -> i > 5) > 1.0
    again_crossed = get_tributary_polygons(_pts(split_outer); holes = [_pts(across)],
                                           hole_supports = [true])
    @test _area(again_crossed, _ -> true) ≈ _area(crossed, _ -> true)
end
