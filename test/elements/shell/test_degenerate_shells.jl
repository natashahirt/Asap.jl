# A triangle with no area makes the shell stiffness undefined. Corners that
# sit on each other become one node. Corners spread along a straight edge
# stay, and that triangle is removed.

using Test
using Asap
using Unitful

_flat_tri(nodes) = Asap.ShellTri3(nodes, 0.2u"m", 30e9u"Pa", 0.2)

@testset "Degenerate shell triangles" begin

    @testset "Corners spread along an edge are not welded" begin
        far1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        far2 = Asap.Node([8.0u"m", 0.0u"m", 0.0u"m"], :free)
        far3 = Asap.Node([24.0u"m", 0.0u"m", 0.0u"m"], :free)
        up = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        over = Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free)
        flat = _flat_tri((far1, far2, far3))
        kept = _flat_tri((far1, up, over))
        @test flat.area < 1e-8
        shells = Asap.ShellTri3[flat, kept]
        panels = [1, 2]
        Asap.collapse_degenerate_shells!(shells; alongside = panels)
        @test length(shells) == 1
        @test panels == [2]
        @test shells[1].nodes === (far1, up, over)
        @test shells[1].area > 1e-8
    end

    @testset "Two corners on the same point become one node" begin
        # q is 0.01 mm from p, on the line to r. The flat triangle welds q
        # onto p, and the triangle that used q follows.
        p = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        q = Asap.Node([1e-5u"m", 0.0u"m", 0.0u"m"], :free)
        r = Asap.Node([8.0u"m", 0.0u"m", 0.0u"m"], :free)
        s = Asap.Node([0.0u"m", 2.0u"m", 0.0u"m"], :free)
        t = Asap.Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)
        flat = _flat_tri((p, q, r))
        uses_q = _flat_tri((q, s, t))
        @test flat.area < 1e-8
        shells = Asap.ShellTri3[flat, uses_q]
        Asap.collapse_degenerate_shells!(shells)
        @test length(shells) == 1
        @test shells[1].nodes[1] === p
        @test shells[1].nodes[2] === s
        @test shells[1].nodes[3] === t
        @test shells[1].area > 1e-6
    end

    @testset "A collinear boundary point is not a shell" begin
        section = Asap.ShellSection(0.2u"m", 30e9u"Pa", 0.2)
        pts = [(0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (0.0, 4.0)]
        corners = Tuple(Asap.Node([pt[1] * u"m", pt[2] * u"m", 0.0u"m"], :free) for pt in pts)
        shells = Asap.Shell(corners, section; target_edge_length = 1.0u"m",
                            edge_support_type = :free, n = 2)
        @test !isempty(shells)
        @test all(el -> el.area > 1e-8, shells)
        @test all(el -> length(unique(objectid, el.nodes)) == 3, shells)
    end

    @testset "alongside must name every shell" begin
        p = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        q = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        r = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        shells = Asap.ShellTri3[_flat_tri((p, q, r))]
        @test_throws ArgumentError Asap.collapse_degenerate_shells!(shells; alongside = Int[])
    end
end
