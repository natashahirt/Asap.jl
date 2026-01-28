using Pkg
# Navigate to Asap project root (test/shell -> test -> Asap)
asap_root = dirname(dirname(@__DIR__))
Pkg.activate(asap_root)
Pkg.instantiate()

include(joinpath(@__DIR__, "test_shell_tri3.jl"))
include(joinpath(@__DIR__, "test_shell_quad4.jl"))
include(joinpath(@__DIR__, "test_mixed_model.jl"))
