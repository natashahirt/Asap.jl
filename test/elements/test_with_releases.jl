using Asap
using Test
using Unitful

@testset "with_releases copies a model with new end releases" begin
    # One 6 m beam between two fully fixed supports under a uniform load w.
    # Fixed-fixed end moments are wL²/12; pinned (:freefree) ends carry none.
    # Both carry wL/2 shear at each end.
    L, w = 6.0, 10_000.0
    n1 = Node([0.0, 0.0, 0.0] .* u"m", :fixed)
    n2 = Node([L, 0.0, 0.0] .* u"m", :fixed)
    section = Section(1.0e-2u"m^2", 30.0u"GPa", 12.0u"GPa",
                      1.0e-4u"m^4", 1.0e-4u"m^4", 1.0e-5u"m^4")
    beam = Element(n1, n2, section, :joist)
    load = LineLoad(beam, [0.0, 0.0, -w] .* u"N/m", :dead)
    case = AbstractLoad[load]
    model = FrameModel([n1, n2], [beam], AbstractLoad[load])

    out = with_releases(model, Dict(1 => :freefree); load_cases = [case])
    copied = out.model

    @test beam isa Element{Asap.FixedFixed}
    @test model.elements[1] === beam
    @test model.loads[1] === load && load.element === beam

    pinned = copied.elements[1]
    @test pinned isa Element{Asap.FreeFree}
    @test pinned !== beam
    @test pinned.id == :joist
    @test pinned.nodeStart === copied.nodes[1] && pinned.nodeEnd === copied.nodes[2]
    @test pinned.nodeStart !== n1
    @test copied.loads[1].element === pinned
    @test copied.loads[1].id == :dead
    @test copied.loads[1] === out.load_cases[1][1]

    solve!(model)
    solve!(copied)
    end_moment(e) = maximum(abs, e.forces[[4, 5, 6, 10, 11, 12]])
    end_shear(e) = maximum(abs, e.forces[[2, 3]])
    @test end_moment(beam) ≈ w * L^2 / 12 rtol = 1e-9
    @test end_moment(pinned) < 1e-6 * w * L^2
    @test end_shear(beam) ≈ w * L / 2 rtol = 1e-9
    @test end_shear(pinned) ≈ w * L / 2 rtol = 1e-9

    @test_throws ArgumentError with_releases(model, Dict(2 => :freefree))
    @test_throws ArgumentError with_releases(model, Dict(1 => :hinged))
end
