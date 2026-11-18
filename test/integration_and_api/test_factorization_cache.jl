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

@testset "A section written without update! is reported stale" begin
    n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2 = Node([3.0u"m", 0.0u"m", 0.0u"m"], :free)
    n3 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :free)
    s = Section(1.0e-2u"m^2", 200.0u"GPa", 77.0u"GPa",
                1.0e-4u"m^4", 1.0e-5u"m^4", 1.0e-6u"m^4")
    e1 = Element(n1, n2, s)
    e2 = Element(n2, n3, s)
    model = FrameModel([n1, n2, n3], [e1, e2],
                       [NodeForce(n3, [0.0u"N", 0.0u"N", -1000.0u"N"])])
    solve!(model)
    @test isempty(stale_frame_elements(model))

    e2.section = Section(s.A, s.E, s.G, 30 * s.Ix, s.Iy, s.J)
    K_before = copy(e2.K)
    @test stale_frame_elements(model) == [2]
    @test e2.K == K_before

    update!(model)
    @test isempty(stale_frame_elements(model))
end

@testset "Springs on stored DOFs keep the symbolic factorization" begin
    # A values-only update drops spring stiffness from S. Putting the springs
    # back on the same DOFs must refactor in place, not from scratch.
    function cantilever(k; E = 100.0u"GPa")
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :free)
        section = Section(
            1.0e-2u"m^2", E, 40.0u"GPa",
            1.0e-4u"m^4", 1.0e-4u"m^4", 1.0e-5u"m^4",
        )
        element = Element(n1, n2, section)
        model = FrameModel([n1, n2], [element],
                           [NodeForce(n2, [0.0u"N", 0.0u"N", -1000.0u"N"])])
        process!(model)
        add_springs!(model, [Spring(n2; kz = k)])
        return model, element, n2
    end

    model, element, n2 = cantilever(1.0e6u"N/m")
    solve!(model)
    factor = model._factorization
    matrix = model.S

    s = element.section
    element.section = Section(s.A, 2 * s.E, s.G, s.Ix, s.Iy, s.J, s.ρ; Ay = s.Ay, Az = s.Az)
    stiffer = Spring(n2; kz = 3.0e6u"N/m")
    empty!(model.springs)
    update!(model; values_only = true)
    add_springs!(model, [stiffer])

    @test model.S === matrix
    @test model._factorization === factor
    @test model._factorization_stale
    solve!(model)
    @test model._factorization === factor

    fresh, _, _ = cantilever(3.0e6u"N/m"; E = 200.0u"GPa")
    solve!(fresh)
    @test model.u ≈ fresh.u rtol = 1e-10
end

@testset "A moved node keeps the symbolic factorization" begin
    # Cantilever along x, tip load downward. Doubling the length scales the
    # tip deflection by 8, and the sparsity pattern does not change.
    n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :free)
    section = Section(
        1.0e-2u"m^2", 100.0u"GPa", 40.0u"GPa",
        1.0e-4u"m^4", 1.0e-4u"m^4", 1.0e-5u"m^4",
    )
    element = Element(n1, n2, section)
    model = FrameModel([n1, n2], [element],
                       [NodeForce(n2, [0.0u"N", 0.0u"N", -1000.0u"N"])])
    solve!(model)
    short = model.u[n2.globalID[3]]
    factor = model._factorization
    matrix = model.S

    n2.position[1] = 4.0u"m"
    update!(model; geometry = true)

    @test model.S === matrix
    @test model._factorization === factor
    @test model._factorization_stale
    solve!(model)
    @test model._factorization === factor
    @test abs(model.u[n2.globalID[3]]) > abs(short) * 2

    n1b = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2b = Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
    fresh = FrameModel([n1b, n2b], [Element(n1b, n2b, section)],
                       [NodeForce(n2b, [0.0u"N", 0.0u"N", -1000.0u"N"])])
    solve!(fresh)
    @test model.u ≈ fresh.u rtol = 1e-8
end
