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
