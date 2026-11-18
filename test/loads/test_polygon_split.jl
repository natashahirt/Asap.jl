# A joist across a trapezoid has to cut two simple pieces. The old walk
# joined the second piece in the wrong order, so it crossed itself and the
# one-way partition counted that area twice.

using Test
using Asap
using Unitful

function _split_beam(x1, y1, x2, y2, id)
    n1 = Asap.Node([x1, y1, 0.0]u"m", :free)
    n2 = Asap.Node([x2, y2, 0.0]u"m", :free)
    sec = Asap.Section(0.005u"m^2", 200e9u"Pa", 77e9u"Pa",
                       1e-4u"m^4", 5e-5u"m^4", 1e-5u"m^4", 7850.0u"kg/m^3")
    return Asap.Element(n1, n2, sec, id)
end

function _shoelace(poly)
    a = 0.0
    n = length(poly)
    for i in 1:n
        x1, y1 = poly[i]
        x2, y2 = poly[mod1(i + 1, n)]
        a += x1 * y2 - x2 * y1
    end
    return abs(a) / 2
end

function _load_area(ld)
    p1 = ld.element.nodeStart.position
    p2 = ld.element.nodeEnd.position
    L = hypot(ustrip(u"m", p2[1] - p1[1]), ustrip(u"m", p2[2] - p1[2]))
    integ = 0.0
    p = ld.positions
    w = ld._widths_m
    for i in 1:(length(p) - 1)
        integ += 0.5 * (w[i] + w[i + 1]) * (p[i + 1] - p[i])
    end
    return L * integ
end

function _joist_share(corners, joist_ends)
    joist = _split_beam(joist_ends... , :joist)
    subs = Asap._split_polygon_by_beams(corners, [joist]; tol = 1e-6)
    tls = Asap._shell_to_tributary_loads_with_interior(
        corners, [joist], [joist], 1.0u"Pa";
        axis = [0.0, 1.0], direction = (0.0, 0.0, -1.0), tol = 1e-4)
    joist_area = sum(_load_area, tls; init = 0.0)
    return (subs = subs, joist_area = joist_area)
end

@testset "a joist cuts a trapezoid into two simple pieces" begin
    narrow = NTuple{2, Float64}[(0.0, 0.0), (8.0, 0.0), (6.0, 4.0), (2.0, 4.0)]
    wide = NTuple{2, Float64}[(0.0, 0.0), (12.0, 0.0), (6.0, 4.0), (2.0, 4.0)]
    cut = _joist_share(narrow, (1.0, 2.0, 7.0, 2.0))
    @test length(cut.subs) == 2
    @test sum(_shoelace, cut.subs) ≈ _shoelace(narrow) atol = 1e-6
    @test all(s -> _shoelace(s) > 1.0, cut.subs)
    @test cut.joist_area > 0.2 * _shoelace(narrow)

    grown = _joist_share(wide, (1.0, 2.0, 9.0, 2.0))
    @test grown.joist_area > cut.joist_area
    @test sum(_shoelace, grown.subs) ≈ _shoelace(wide) atol = 1e-6
end
