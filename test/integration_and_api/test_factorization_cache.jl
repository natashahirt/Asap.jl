using Asap
using Test
using Unitful

@testset "Factorization cache invalidation" begin
    # A values-only stiffness mutation must mark the cached factor stale and
    # refactor it before the next solve, without rebuilding the sparsity pattern.
    n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :xfree)
    section = Section(
        1.0e-2u"m^2", 100.0u"GPa", 40.0u"GPa",
        1.0e-4u"m^4", 1.0e-4u"m^4", 1.0e-5u"m^4",
    )
    element = Element(n1, n2, section)
    model = FrameModel(
        [n1, n2],
        [element],
        [NodeForce(n2, [1000.0u"N", 0.0u"N", 0.0u"N"])],
    )

    solve!(model)
    initial_displacement = model.u[n2.globalID[1]]
    initial_factor = model._factorization
    initial_matrix = model.S

    element.section = Section(
        section.A, 2 * section.E, section.G,
        section.Ix, section.Iy, section.J, section.ρ;
        Ay=section.Ay, Az=section.Az,
    )
    update!(model; values_only=true)

    @test model.S === initial_matrix
    @test model._factorization === initial_factor
    @test model._factorization_stale

    solve!(model)

    @test !model._factorization_stale
    @test model._factorization === initial_factor
    @test model.u[n2.globalID[1]] ≈ initial_displacement / 2
end
