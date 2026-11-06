#=
Constrained Delaunay Interior Support Lines
===========================================

Tests for `Shell(...; interior_constrained_lines=...)` (Phase A3).

Unlike `interior_supports` (sampled as points + given boundary conditions),
`interior_constrained_lines` makes each consecutive-vertex pair of a polyline a
*constrained edge* of the triangulation — the mesh is guaranteed to have an
element edge aligned with each segment. This is what lets a slab FEA mesh
respect joist / framing lines exactly rather than by point sampling.

Edge presence is checked against the RETURNED shell mesh (not the intermediate
triangulation): in a triangle every pair of its three vertices is a mesh edge,
so if the two polyline vertices are both corners of some `ShellTri3`, the
segment survives as an edge.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

# Register Asap's custom units when run standalone (harmless if already registered)
try
    Unitful.register(Asap)
catch
end

"""Index of the shell-mesh node nearest `target=(x,y)`; nothing if none within tol."""
function _node_at(shells, target::Tuple{Float64, Float64}; tol::Float64 = 1e-6)
    for s in shells
        for nd in s.nodes
            x = ustrip(u"m", nd.position[1])
            y = ustrip(u"m", nd.position[2])
            if hypot(x - target[1], y - target[2]) < tol
                return nd
            end
        end
    end
    return nothing
end

"""
True if `a` and `b` are both vertices of a common triangle in `shells` — i.e.
the segment `a–b` is an edge of the returned mesh.
"""
function _mesh_edge_present(shells, a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64};
                            tol::Float64 = 1e-6)
    for s in shells
        found_a = false
        found_b = false
        for nd in s.nodes
            x = ustrip(u"m", nd.position[1])
            y = ustrip(u"m", nd.position[2])
            hypot(x - a[1], y - a[2]) < tol && (found_a = true)
            hypot(x - b[1], y - b[2]) < tol && (found_b = true)
        end
        found_a && found_b && return true
    end
    return false
end

@testset "Constrained Delaunay interior lines" begin

    sec = ShellSection(0.15u"m", 30u"GPa", 0.2)

    @testset "Rectangular slab — one interior support line" begin
        # 6m x 4m slab. A bent (non-axis-aligned) polyline so we exercise two
        # consecutive segments and avoid accidental grid-line coincidence.
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)

        line = [(1.0, 1.0), (3.0, 2.5), (5.0, 1.0)]  # NTuple{2,Float64} polyline

        shells = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            interior_constrained_lines = [line],
        )

        @test length(shells) > 0

        # Every polyline vertex must be a mesh node …
        for v in line
            @test _node_at(shells, v) !== nothing
        end
        # … and every consecutive vertex pair must be a mesh edge.
        for i in 1:(length(line) - 1)
            @test _mesh_edge_present(shells, line[i], line[i + 1])
        end
    end

    @testset "L-shaped slab — multiple support lines" begin
        # Concave L: x∈[0,6],y∈[0,3] plus x∈[0,3],y∈[3,6].
        c = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned),
            Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned),
            Node([6.0u"m", 3.0u"m", 0.0u"m"], :pinned),
            Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned),
            Node([3.0u"m", 6.0u"m", 0.0u"m"], :pinned),
            Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned),
        ]

        line_a = [(1.0, 1.0), (3.0, 1.5), (5.0, 1.0)]  # in the wide bottom leg
        line_b = [(1.0, 4.0), (2.0, 4.5)]              # in the tall left leg

        shells = Shell(c, sec;
            target_edge_length = 1.0u"m",
            interior_constrained_lines = [line_a, line_b],
        )

        @test length(shells) > 0

        for line in (line_a, line_b)
            for v in line
                @test _node_at(shells, v) !== nothing
            end
            for i in 1:(length(line) - 1)
                @test _mesh_edge_present(shells, line[i], line[i + 1])
            end
        end
    end

    @testset "Support line through a column patch" begin
        # 6m x 4m slab with a stiff column patch centered at (3, 2), and a
        # constrained line running horizontally through the patch interior.
        # The patch vertices sit at y = 1.7 / 2.3, so none is collinear with the
        # y = 2.0 line — the constrained edge must survive and the mesh must build.
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)

        stiff = ShellSection(0.5u"m", 200e9u"Pa", 0.3; name = :column)
        patch = ShellPatch(3.0, 2.0, 0.6, 0.6, stiff; id = :col_patch)

        line = [(1.0, 2.0), (5.0, 2.0)]

        shells = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            interior_patches = [patch],
            interior_constrained_lines = [line],
        )

        @test length(shells) > 0
        # Constrained edge preserved despite passing through the patch region.
        @test _mesh_edge_present(shells, line[1], line[2]) ||
              # If a patch/grid vertex lands on the interior, the edge may split;
              # then both sub-edges through the shared midpoint must exist.
              (_node_at(shells, (3.0, 2.0)) !== nothing &&
               _mesh_edge_present(shells, line[1], (3.0, 2.0)) &&
               _mesh_edge_present(shells, (3.0, 2.0), line[2]))
        # Patch section was applied to at least one element.
        @test any(s -> s.id == :col_patch, shells)
    end

    @testset "Backward compatibility — empty lines" begin
        # With no constrained lines the result must match the pre-existing path.
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)

        shells_default = Shell((n1, n2, n3, n4), sec; target_edge_length = 1.0u"m")
        shells_empty = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            interior_constrained_lines = Vector{NTuple{2, Float64}}[],
        )

        @test length(shells_default) == length(shells_empty)
        @test length(shells_default) > 0
    end
end
