using Asap
using SparseArrays
using Test
using Unitful

const _PP_SECTION = Section(
    0.01u"m^2", 200.0u"GPa", 80.0u"GPa",
    1.0e-4u"m^4", 1.0e-4u"m^4", 1.0e-5u"m^4",
)

function _parametric_fixture(
    model_type;
    host_release::Symbol=:fixedfixed,
    loads_builder=(_host_a, _bridge, _nodes) -> AbstractLoad[],
)
    nodes = [
        Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
        Node([10.0u"m", 0.0u"m", 0.0u"m"], :free),
        Node([0.0u"m", 4.0u"m", 0.0u"m"], :fixed),
        Node([10.0u"m", 4.0u"m", 0.0u"m"], :free),
    ]
    host_a = Element(
        nodes[1], nodes[2], _PP_SECTION, :host_a;
        release=host_release,
        eccentricity=(0.1u"m", 0.2u"m", 0.03u"m", 0.04u"m"),
    )
    host_b = Element(nodes[3], nodes[4], _PP_SECTION, :host_b)
    host_a.Ψ = 0.17
    bridge = ParametricFrameElement(
        host_a, 0.25, host_b, 0.75, _PP_SECTION, :bridge;
        release=:fixedfree,
    )
    bridge.Ψ = 0.31
    elements = Asap.FrameElement[host_a, host_b, bridge]
    loads = AbstractLoad[loads_builder(host_a, bridge, nodes)...]
    model = if model_type === FrameModel
        FrameModel(nodes, elements, loads)
    else
        shell = ShellTri3(
            (nodes[1], nodes[3], nodes[4]), 0.1u"m", 25.0u"GPa", 0.2;
            id=:retained_shell,
        )
        Model(nodes, elements, ShellElement[shell], loads)
    end
    return (; model, nodes, host_a, host_b, bridge)
end

_frames(model::FrameModel) = model.elements
_frames(model::Model) = model.frame_elements

@testset "Parametric frame processing" begin
    for model_type in (FrameModel, Model)
        @testset "$(nameof(model_type)) topology and load remapping" begin
            fixture = _parametric_fixture(
                model_type;
                loads_builder=(host, bridge, nodes) -> AbstractLoad[
                    NodeForce(nodes[1], [1.0u"N", 2.0u"N", 3.0u"N"], :node),
                    LineLoad(host, [0.0u"N/m", -2.0u"N/m", 0.0u"N/m"], :line),
                    PointLoad(host, 0.10, [0.0u"N", -10.0u"N", 0.0u"N"], :before),
                    PointLoad(host, 0.25, [0.0u"N", -20.0u"N", 0.0u"N"], :at),
                    PointLoad(host, 0.60, [0.0u"N", -30.0u"N", 0.0u"N"], :after),
                    LineLoad(bridge, [0.0u"N/m", 0.0u"N/m", -3.0u"N/m"], :bridge_line),
                    PointLoad(bridge, 0.40, [0.0u"N", 0.0u"N", -4.0u"N"], :bridge_point),
                    NodeMoment(nodes[1], [1.0u"N*m", 2.0u"N*m", 3.0u"N*m"], :moment),
                ],
            )

            retained_shell = model_type === Model ? only(fixture.model.shell_elements) : nothing
            process!(fixture.model)
            frames = _frames(fixture.model)

            @test length(fixture.model.nodes) == 6
            @test length(frames) == 5
            @test !any(element -> element isa ParametricFrameElement, frames)
            @test fixture.model.nNodes == 6
            @test size(fixture.model.S) == (36, 36)
            @test fixture.model.nDOFs == 36
            @test sort([node.nodeID for node in fixture.model.nodes]) == collect(1:6)
            @test sort([element.elementID for element in frames]) == collect(1:5)
            if model_type === Model
                @test only(fixture.model.shell_elements) === retained_shell
                @test fixture.model.nShellElements == 1
            end

            host_a_fragments = [element for element in frames if element.id === :host_a]
            host_b_fragments = [element for element in frames if element.id === :host_b]
            resolved_bridge = only([element for element in frames if element.id === :bridge])
            @test sort(ustrip.(u"m", getproperty.(host_a_fragments, :length))) ≈ [2.5, 7.5]
            @test sort(ustrip.(u"m", getproperty.(host_b_fragments, :length))) ≈ [2.5, 7.5]
            @test resolved_bridge isa Element{FixedFree}
            @test resolved_bridge.section === _PP_SECTION
            @test resolved_bridge.Ψ == 0.31
            @test fixture.host_a.nodeStart === fixture.nodes[1]
            @test fixture.host_a.nodeEnd === fixture.nodes[2]
            @test fixture.host_a.elementID == 0

            ordered_a = sort(host_a_fragments; by=e -> to_meters(e.nodeStart.position[1]))
            @test all(element -> element.section === _PP_SECTION, ordered_a)
            @test all(element -> element.Ψ == 0.17, ordered_a)
            @test ordered_a[1].eccentricity == (0.1, 0.0, 0.03, 0.04)
            @test ordered_a[2].eccentricity == (0.0, 0.2, 0.03, 0.04)

            line_parts = [load for load in fixture.model.loads if load.id === :line]
            @test length(line_parts) == 2
            @test sum(to_meters(load.element.length) * ustrip(u"N/m", load.value[2])
                      for load in line_parts) ≈ -20.0
            @test only([load for load in fixture.model.loads if load.id === :at]) isa NodeForce
            @test only([load for load in fixture.model.loads if load.id === :before]).position ≈ 0.4
            @test only([load for load in fixture.model.loads if load.id === :after]).position ≈ 7 / 15
            @test only([load for load in fixture.model.loads if load.id === :bridge_line]).element ===
                  resolved_bridge
            @test only([load for load in fixture.model.loads if load.id === :bridge_point]).element ===
                  resolved_bridge
            @test count(load -> load.id === :node, fixture.model.loads) == 1
            @test only([load for load in fixture.model.loads if load.id === :moment]) isa NodeMoment
            @test ustrip.(u"N", only([
                load for load in fixture.model.loads if load.id === :at
            ]).value) == [0.0, -20.0, 0.0]
        end

        @testset "$(nameof(model_type)) repeat processing and springs" begin
            fixture = _parametric_fixture(model_type)
            process!(fixture.model)
            spring = Spring(fixture.nodes[2]; kx=1.0e6u"N/m", id=:endpoint_spring)
            add_springs!(fixture.model, spring)
            stiffness_with_spring = copy(fixture.model.S)
            add_springs!(fixture.model, spring)
            @test fixture.model.S == stiffness_with_spring
            signature_1 = (
                length(fixture.model.nodes),
                length(_frames(fixture.model)),
                copy(fixture.model.P),
                copy(fixture.model.Pf),
                copy(fixture.model.S),
            )

            process!(fixture.model)
            signature_2 = (
                length(fixture.model.nodes),
                length(_frames(fixture.model)),
                copy(fixture.model.P),
                copy(fixture.model.Pf),
                copy(fixture.model.S),
            )
            @test signature_2 == signature_1
            @test fixture.model.S[fixture.nodes[2].globalID[1], fixture.nodes[2].globalID[1]] >= 1.0e6
            @test length(fixture.model.springs) == 1
        end
    end

    @testset "coincident landings are tolerance-merged" begin
        fixture = _parametric_fixture(FrameModel)
        second = ParametricFrameElement(
            fixture.host_a, 0.25 + fixture.model.tol / 2,
            fixture.host_b, 0.75 + fixture.model.tol / 2,
            _PP_SECTION, :bridge_2,
        )
        push!(fixture.model.elements, second)
        process!(fixture.model)
        @test count(node -> node.id === :parametric, fixture.model.nodes) == 2
        @test count(element -> element.id in (:bridge, :bridge_2), fixture.model.elements) == 2
    end

    @testset "springs on generated nodes fail before replacement" begin
        fixture = _parametric_fixture(FrameModel)
        process!(fixture.model)
        generated = only([
            node for node in fixture.model.nodes
            if node.id === :parametric && isapprox(to_meters(node.position[2]), 0.0)
        ])
        add_springs!(fixture.model, Spring(generated; kz=1.0e6u"N/m"))
        nodes_before = copy(fixture.model.nodes)
        elements_before = copy(fixture.model.elements)
        @test_throws ArgumentError process!(fixture.model)
        @test fixture.model.nodes == nodes_before
        @test fixture.model.elements == elements_before
        @test !fixture.model.processed
    end

    @testset "release preservation" begin
        expected = Dict(
            :fixedfixed => (FixedFixed, FixedFixed),
            :fixedfree => (FixedFixed, FixedFree),
            :freefixed => (FreeFixed, FixedFixed),
            :freefree => (FreeFixed, FixedFree),
        )
        for (release, types) in expected
            fixture = _parametric_fixture(FrameModel; host_release=release)
            process!(fixture.model)
            fragments = sort(
                [element for element in fixture.model.elements if element.id === :host_a];
                by=e -> to_meters(e.nodeStart.position[1]),
            )
            @test fragments[1] isa Element{types[1]}
            @test fragments[2] isa Element{types[2]}
        end
        fixture = _parametric_fixture(FrameModel; host_release=:joist)
        @test_throws ArgumentError process!(fixture.model)
        @test !fixture.model.processed
    end

    @testset "unsupported loads fail atomically" begin
        fixture = _parametric_fixture(
            FrameModel;
            loads_builder=(host, _bridge, _nodes) -> AbstractLoad[
                TributaryLoad(
                    host, [0.0, 1.0], [1.0u"m", 1.0u"m"], 2.0u"kPa"),
            ],
        )
        original_elements = copy(fixture.model.elements)
        original_nodes = copy(fixture.model.nodes)
        @test_throws ArgumentError process!(fixture.model)
        @test fixture.model.elements == original_elements
        @test fixture.model.nodes == original_nodes
        @test !fixture.model.processed
    end

    @testset "unit-equivalent topology" begin
        n1 = Node([0.0u"ft", 0.0u"ft", 0.0u"ft"], :fixed)
        n2 = Node([10.0u"ft", 0.0u"ft", 0.0u"ft"], :free)
        n3 = Node([0.0u"ft", 4.0u"ft", 0.0u"ft"], :fixed)
        n4 = Node([10.0u"ft", 4.0u"ft", 0.0u"ft"], :free)
        a = Element(n1, n2, _PP_SECTION, :a)
        b = Element(n3, n4, _PP_SECTION, :b)
        p = ParametricFrameElement(a, 0.25, b, 0.75, _PP_SECTION)
        model = FrameModel(
            [n1, n2, n3, n4],
            Asap.FrameElement[a, b, p],
            AbstractLoad[LineLoad(a, [0.0u"kip/ft", -1.0u"kip/ft", 0.0u"kip/ft"])],
        )
        process!(model)
        fracture = only([
            node for node in model.nodes
            if node.id === :parametric && isapprox(to_meters(node.position[2]), 0.0)
        ])
        @test to_meters(fracture.position[1]) ≈ 0.762
    end
end

