using Asap
using StructuralBase: StructuralUnits
using Test
using LinearAlgebra
using Unitful

Unitful.register(StructuralUnits)

@testset "Asap Tests" begin
    include("core.jl")
    include("tributary_load/test_tributary_load.jl")
    include("shell/test_mixed_model.jl")
    include("shell/test_shell_quad4.jl")
    include("shell/test_shell_tri3.jl")
end
