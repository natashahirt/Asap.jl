# Tributary partitions: every square meter of floor lands on exactly one edge,
# and each edge's load profile carries exactly its area.

using Test
using Asap
using Meshes: Point
using Unitful: @u_str

_pts(xy) = [Point(x * u"m", y * u"m") for (x, y) in xy]
_len(xy, i) = hypot(xy[mod1(i + 1, length(xy))][1] - xy[i][1], xy[mod1(i + 1, length(xy))][2] - xy[i][2])
_areas(tribs) = Dict(t.local_edge_idx => t.area for t in sort(tribs; by = t -> t.local_edge_idx))
_carried(t, xy) = Asap.profile_area(t, _len(xy, t.local_edge_idx))

@testset "Tributary balance" begin
    rect = [(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)]

    @testset "two-way rectangle: 45° trapezoids and triangles" begin
        tribs = get_tributary_polygons_isotropic(_pts(rect))
        a = _areas(tribs)
        @test a[1] ≈ 8.0 atol = 1e-6
        @test a[2] ≈ 4.0 atol = 1e-6
        @test a[3] ≈ 8.0 atol = 1e-6
        @test a[4] ≈ 4.0 atol = 1e-6
        for t in tribs
            @test _carried(t, rect) ≈ t.area atol = 1e-9
            @test t.profile_s[1] == 0.0 && t.profile_s[end] == 1.0
        end
        bottom = tribs[findfirst(t -> t.local_edge_idx == 1, tribs)]
        @test maximum(bottom.profile_w) ≈ 2.0 atol = 1e-6
    end

    @testset "two-way collinear chain that wraps past vertex 1" begin
        # Bottom side is edge 5 (0,0)→(2,0) then edge 1 (2,0)→(6,0).
        xy = [(2.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0), (0.0, 0.0)]
        tribs = get_tributary_polygons_isotropic(_pts(xy))
        a = _areas(tribs)
        @test a[5] ≈ 2.0 atol = 1e-6   # under y = x for x ∈ [0, 2]
        @test a[1] ≈ 6.0 atol = 1e-6
        for t in tribs
            @test all(s -> -1e-9 <= s <= 1 + 1e-9, t.s)
            @test _carried(t, xy) ≈ t.area atol = 1e-9
        end
    end

    @testset "a free edge carries no one-way load" begin
        # 6 m by 4 m, spanning in x. Equal supports split the 24 m². Removing
        # the left edge's support moves all of it onto the right edge, and
        # removing both supports is refused.
        tribs = get_tributary_polygons(_pts(rect); axis = [1.0, 0.0])
        full = get_tributary_polygons(_pts(rect); axis = [1.0, 0.0], weights = [1.0, 1.0, 1.0, 0.0])
        @test _areas(tribs)[2] ≈ 12.0 atol = 1e-6
        @test _areas(full)[2] ≈ 24.0 atol = 1e-6
        @test _areas(full)[4] ≈ 0.0 atol = 1e-6
        @test sum(t.area for t in full) ≈ 24.0 atol = 1e-6
        @test_throws ArgumentError get_tributary_polygons(
            _pts(rect); axis = [1.0, 0.0], weights = [0.0, 0.0, 0.0, 0.0])
    end

    @testset "one-way onto skewed edges keeps the whole strip" begin
        par = [(0.0, 0.0), (5.0, 0.0), (6.5, 3.0), (1.5, 3.0)]
        tribs = get_tributary_polygons(_pts(par); axis = [1.0, 0.0])
        a = _areas(tribs)
        @test a[2] ≈ 7.5 atol = 1e-9
        @test a[4] ≈ 7.5 atol = 1e-9
        @test a[1] ≈ 0.0 atol = 1e-12
        right = tribs[findfirst(t -> t.local_edge_idx == 2, tribs)]
        @test _carried(right, par) ≈ 7.5 atol = 1e-9
        # 2.5 m of strip meets an edge at sin θ = 3 / √11.25
        @test all(w -> isapprox(w, 2.5 * 3 / sqrt(11.25); atol = 1e-9), right.profile_w)

        tribs_y = get_tributary_polygons(_pts(par); axis = [0.0, 1.0])
        @test sum(t.area for t in tribs_y) ≈ 15.0 atol = 1e-9
        for t in tribs_y
            @test _carried(t, par) ≈ t.area atol = 1e-9
        end
    end

    slab = [(0.0, 0.0), (10.0, 0.0), (10.0, 4.0), (0.0, 4.0)]
    hole = [(4.0, 1.0), (6.0, 1.0), (6.0, 3.0), (4.0, 3.0)]

    @testset "one-way: an opening edge carries nothing" begin
        tribs = get_tributary_polygons(_pts(slab); axis = [0.0, 1.0], holes = [hole])
        a = _areas(tribs)
        @test a[1] ≈ 18.0 atol = 1e-9
        @test a[3] ≈ 18.0 atol = 1e-9
        @test a[2] ≈ 0.0 atol = 1e-12
        @test a[4] ≈ 0.0 atol = 1e-12
        bottom = tribs[findfirst(t -> t.local_edge_idx == 1, tribs)]
        @test _carried(bottom, slab) ≈ 18.0 atol = 1e-9
        @test minimum(bottom.profile_w) ≈ 1.0 atol = 1e-9   # behind the opening
        @test maximum(bottom.profile_w) ≈ 2.0 atol = 1e-9
    end

    @testset "one-way: opening across the boundary sends its strip to the far edge" begin
        notch = [(4.0, -1.0), (6.0, -1.0), (6.0, 1.0), (4.0, 1.0)]
        tribs = get_tributary_polygons(_pts(slab); axis = [0.0, 1.0], holes = [notch])
        a = _areas(tribs)
        @test a[1] ≈ 16.0 atol = 1e-9
        @test a[3] ≈ 22.0 atol = 1e-9
    end

    @testset "one-way: a strip between two openings has no support" begin
        low = [(4.0, 0.5), (6.0, 0.5), (6.0, 1.5), (4.0, 1.5)]
        high = [(4.0, 2.5), (6.0, 2.5), (6.0, 3.5), (4.0, 3.5)]
        @test_throws ArgumentError get_tributary_polygons(_pts(slab); axis = [0.0, 1.0], holes = [low, high])
    end

    @testset "two-way: each region loses the part inside the opening" begin
        tribs = get_tributary_polygons_isotropic(_pts(rect); holes = [[(2.0, 1.0), (4.0, 1.0), (4.0, 3.0), (2.0, 3.0)]])
        a = _areas(tribs)
        @test a[1] ≈ 6.0 atol = 1e-6
        @test a[2] ≈ 4.0 atol = 1e-6
        @test a[3] ≈ 6.0 atol = 1e-6
        @test a[4] ≈ 4.0 atol = 1e-6
        for t in tribs
            @test _carried(t, rect) ≈ t.area atol = 1e-9
        end
        @test sum(t.fraction for t in tribs) ≈ 1.0 atol = 1e-6
    end

    @testset "a face its opening fills carries nothing" begin
        for axis in (nothing, [1.0, 0.0])
            tribs = get_tributary_polygons(_pts(rect); axis = axis, holes = [rect])
            @test sum(t.area for t in tribs) ≈ 0.0 atol = 1e-9
        end
    end

    @testset "two-way re-entrant corner: the wedge past an edge end folds back onto it" begin
        ell = [(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (6.0, 6.0), (6.0, 2.0), (0.0, 2.0)]
        tribs = get_tributary_polygons_isotropic(_pts(ell))
        @test sum(t.area for t in tribs) ≈ 36.0 atol = 1e-6   # 10×2 + 4×4
        for t in tribs
            @test _carried(t, ell) ≈ t.area atol = 1e-6
        end
        # Edge 4 runs (6,6) → (6,2) and ends at the re-entrant corner.
        corner = tribs[findfirst(t -> t.local_edge_idx == 4, tribs)]
        @test any(s -> s > 1 + 1e-9, corner.s)
        @test corner.profile_w[end - 1] > 0.0
    end

    @testset "two-way split events: a reflex vertex cuts the wavefront in two" begin
        # The slot floor's front meets the bottom's front at t = 0.5, between
        # x = 3 - t and 5 + t: ∫₀^½ (2 + 2t) dt = 1.25. Each wall translates
        # 0.5 m, then takes a 0.25 m² end triangle in its 2 × 1 room.
        slot = [(0.0, 0.0), (8.0, 0.0), (8.0, 2.0), (5.0, 2.0), (5.0, 1.0), (3.0, 1.0),
                (3.0, 2.0), (0.0, 2.0)]
        a = _areas(get_tributary_polygons_isotropic(_pts(slot)))
        @test a[5] ≈ 1.25 atol = 1e-6
        @test a[4] ≈ 0.75 atol = 1e-6
        @test a[6] ≈ 0.75 atol = 1e-6
        @test sum(values(a)) ≈ 14.0 atol = 1e-6

        # Two apexes meet head-on (a vertex event); the result keeps both symmetries.
        notches = [(0.0, 0.0), (5.0, 0.0), (6.0, 1.5), (7.0, 0.0), (12.0, 0.0), (12.0, 4.0),
                   (7.0, 4.0), (6.0, 2.5), (5.0, 4.0), (0.0, 4.0)]
        tribs = get_tributary_polygons_isotropic(_pts(notches))
        a = _areas(tribs)
        @test sum(values(a)) ≈ 45.0 atol = 1e-6
        @test a[1] ≈ a[4] ≈ a[6] ≈ a[9]
        @test a[2] ≈ a[3] ≈ a[7] ≈ a[8]
        @test a[5] ≈ 4.0 atol = 1e-6
        for t in tribs
            @test _carried(t, notches) ≈ t.area atol = 1e-6
        end

        # A rib bay cut by an angled grid at an opening corner (4, 2).
        bay = [(4.0, 2.0), (4.0, 2.403), (3.6662, 2.2103), (4.1162, 1.4309),
               (4.8956, 1.8809), (4.8268, 2.0)]
        tribs = get_tributary_polygons_isotropic(_pts(bay))
        @test sum(t.area for t in tribs) ≈ Asap._net_area(bay, Vector{NTuple{2,Float64}}[]) atol = 1e-9
    end

    @testset "two-way: a region past a short edge folds until it lies on the edge" begin
        # Edge 1 is 0.15 m long and starts at the re-entrant corner (2, 2);
        # its region runs past that corner by more than its own length.
        bay = [(2.0, 2.0), (2.0, 2.1485), (1.8457, 2.5724), (1.0, 2.2646), (1.3078, 1.4189),
               (2.1535, 1.7267), (2.0541, 2.0)]
        tribs = get_tributary_polygons_isotropic(_pts(bay))
        short = tribs[findfirst(t -> t.local_edge_idx == 1, tribs)]
        @test any(s -> s < -1.0, short.s)
        for t in tribs
            @test _carried(t, bay) ≈ t.area atol = 1e-9
        end
    end

    @testset "profile sum keeps steps" begin
        p = Asap._sum_profiles([([0.0, 0.5, 0.5, 1.0], [1.0, 1.0, 3.0, 3.0]),
                                ([0.0, 1.0], [1.0, 1.0])])
        @test Asap._trapezoid(p...) ≈ 3.0 atol = 1e-12
    end

    @testset "hand-built outline gets a sliced profile" begin
        t = TributaryPolygon(1, [0.0, 1.0, 0.75, 0.25], [0.0, 0.0, 1.0, 1.0], 3.0, 0.5)
        @test Asap.profile_area(t, 4.0) ≈ 3.0 atol = 1e-12
    end
end
