using Asap
using StructuralBase: StructuralUnits
using Test
using LinearAlgebra
using Unitful

Unitful.register(StructuralUnits)

@testset "Asap Tests" begin
    # Core functionality
    include("core.jl")
    
    # =========================================================================
    # Elements
    # =========================================================================
    @testset "Elements" begin
        include("elements/test_beam_features.jl")
        include("elements/test_spring.jl")
        
        # Shell elements
        include("elements/shell/test_shell_tri3.jl")
        include("elements/shell/test_composite_shell.jl")
        include("elements/shell/test_benchmarks.jl")
        include("elements/shell/test_plate_analytical.jl")
        include("elements/shell/test_finetools_validation.jl")
        include("elements/shell/test_scordelis_lo.jl")
        include("elements/shell/test_meshing.jl")
        include("elements/shell/test_interior_supports.jl")
        include("elements/shell/test_diaphragm.jl")
    end
    
    # =========================================================================
    # Loads
    # =========================================================================
    @testset "Loads" begin
        include("loads/test_tributary_load.jl")
        include("loads/test_area_load.jl")
    end
    
    # =========================================================================
    # Tributary Geometry
    # =========================================================================
    @testset "Tributary Geometry" begin
        include("tributary/test_spans.jl")
        include("tributary/test_voronoi_tributaries.jl")
        # Note: test_tributaries.jl and test_cell_depths.jl are debug/visualization scripts
        # that require GLMakie and should be run manually
    end
    
    # =========================================================================
    # Integration & API
    # =========================================================================
    @testset "Integration & API" begin
        include("integration_and_api/test_mixed_model.jl")
        include("integration_and_api/test_shell_beam_interactions.jl")
        include("integration_and_api/test_shell_internal_forces.jl")
        include("integration_and_api/test_new_api.jl")
        include("integration_and_api/test_load_approaches.jl")
    end
    
    # =========================================================================
    # Nonlinear Analysis
    # =========================================================================
    @testset "Nonlinear Analysis" begin
        include("nonlinear/test_buckling.jl")
        include("nonlinear/test_pdelta.jl")
        include("nonlinear/test_nonlinear_statics.jl")
    end
end
