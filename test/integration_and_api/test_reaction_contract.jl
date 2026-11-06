using Asap
using Test
using Unitful

function _restrained_frame_model(model_type)
    n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
    section = Section(
        0.01u"m^2", 200.0u"GPa", 80.0u"GPa",
        1.0e-4u"m^4", 1.0e-4u"m^4", 1.0e-5u"m^4",
    )
    element = Element(n1, n2, section)
    load = NodeForce(n1, [10.0u"N", -20.0u"N", 30.0u"N"])
    return model_type([n1, n2], [element], [load]), n1
end

@testset "Restrained-DOF nodal-load reactions" begin
    @testset "FrameModel" begin
        model, loaded_support = _restrained_frame_model(FrameModel)
        solve!(model)
        dofs = loaded_support.globalID[1:3]
        @test model.reactions[dofs] ≈ [-10.0, 20.0, -30.0]
    end

    @testset "mixed Model" begin
        model, loaded_support = _restrained_frame_model(Model)
        solve!(model)
        dofs = loaded_support.globalID[1:3]
        @test model.reactions[dofs] ≈ [-10.0, 20.0, -30.0]
    end

    @testset "TrussModel" begin
        n1 = TrussNode([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = TrussNode([1.0u"m", 0.0u"m", 0.0u"m"], :xfree)
        element = TrussElement(n1, n2, TrussSection(0.01u"m^2", 200.0u"GPa"))
        load = NodeForce(n1, [10.0u"N", -20.0u"N", 30.0u"N"])
        model = TrussModel([n1, n2], [element], [load])

        solve!(model)

        @test model.reactions[n1.globalID] ≈ [-10.0, 20.0, -30.0]
    end
end
