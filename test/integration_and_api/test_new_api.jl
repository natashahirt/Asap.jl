#=
Tests for New API Features
===========================

Tests for ShellMaterial, modal!/modal, and other new API additions.
These tests verify both new and old approaches work correctly.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "New API Features" begin
    
    @testset "ShellMaterial" begin
        @testset "Construction" begin
            # Basic construction
            concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
            @test concrete.E ≈ 30e9
            @test concrete.ν ≈ 0.2
            @test concrete.ρ ≈ 2400.0
            @test concrete.name == :concrete
            
            # Default density
            mat = ShellMaterial(E=200u"GPa", ν=0.3)
            @test mat.ρ ≈ 0.0
            @test mat.name == :material
            
            # Unit conversion
            mat_ksi = ShellMaterial(E=4000u"ksi", ν=0.2, ρ=150u"lb/ft^3")
            @test mat_ksi.E ≈ ustrip(u"Pa", 4000u"ksi") rtol=1e-6
        end
        
        @testset "Inline Material Creation" begin
            # Presets removed - materials created inline
            concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
            @test concrete.name == :concrete
            @test concrete.E ≈ 30e9
            @test concrete.ν ≈ 0.2
            @test concrete.ρ ≈ 2400.0
            
            steel = ShellMaterial(E=200u"GPa", ν=0.3, ρ=7850u"kg/m^3", name=:steel)
            @test steel.name == :steel
            @test steel.E ≈ 200e9
            @test steel.ν ≈ 0.3
            @test steel.ρ ≈ 7850.0
        end
        
        @testset "Invalid Poisson Ratio" begin
            @test_throws AssertionError ShellMaterial(E=30u"GPa", ν=0.5)
            @test_throws AssertionError ShellMaterial(E=30u"GPa", ν=0.0)
            @test_throws AssertionError ShellMaterial(E=30u"GPa", ν=-0.1)
        end
    end
    
    @testset "ShellTri3 with ShellMaterial" begin
        # Create test nodes
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)
        
        @testset "Old API still works" begin
            # Original constructor: (nodes, thickness, E, ν; ρ)
            shell_old = ShellTri3((n1, n2, n3), 0.2u"m", 30e9u"Pa", 0.2; ρ=2400.0)
            @test shell_old.thickness ≈ 0.2
            @test shell_old.E ≈ 30e9
            @test shell_old.ν ≈ 0.2
            @test shell_old.ρ ≈ 2400.0
        end
        
        @testset "New API with ShellMaterial" begin
            concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
            shell_new = ShellTri3((n1, n2, n3), 0.2u"m", concrete)
            
            @test shell_new.thickness ≈ 0.2
            @test shell_new.E ≈ 30e9
            @test shell_new.ν ≈ 0.2
            @test shell_new.ρ ≈ 2400.0
        end
        
        @testset "Both approaches give same stiffness" begin
            concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3")
            
            # Reset nodes for fresh shells
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3 = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)
            
            shell_old = ShellTri3((n1, n2, n3), 0.2u"m", 30e9u"Pa", 0.2; ρ=2400.0)
            
            n1b = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2b = Node([1.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3b = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)
            
            shell_new = ShellTri3((n1b, n2b, n3b), 0.2u"m", concrete)
            
            # Process both to populate stiffness matrices
            model_old = ShellModel([n1, n2, n3], [shell_old], AbstractLoad[])
            model_new = ShellModel([n1b, n2b, n3b], [shell_new], AbstractLoad[])
            
            process!(model_old)
            process!(model_new)
            
            # Stiffness matrices should be identical
            @test shell_old.K ≈ shell_new.K
        end
        
        @testset "Iteration with ShellMaterial" begin
            # This tests the main use case: iterating over materials
            materials = [
                ShellMaterial(E=30u"GPa", ν=0.2, name=:concrete),
                ShellMaterial(E=200u"GPa", ν=0.3, name=:steel),
            ]
            
            results = Dict{Symbol, Float64}()
            
            for mat in materials
                n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
                n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :fixed)
                n3 = Node([0.5u"m", 0.5u"m", 0.0u"m"], :free)
                
                shell = ShellTri3((n1, n2, n3), 0.1u"m", mat)
                load = NodeForce(n3, [0.0u"N", 0.0u"N", -1000.0u"N"])
                
                model = ShellModel([n1, n2, n3], [shell], [load])
                solve!(model)
                
                # Store deflection
                results[mat.name] = abs(ustrip(u"m", n3.displacement[3]))
            end
            
            # Steel should deflect less (stiffer)
            @test results[:steel] < results[:concrete]
        end
    end
    
    @testset "Modal API: modal! and modal" begin
        # Create a simple cantilever beam for modal testing
        L = 2.0u"m"
        E = 200u"GPa"
        ρ = 7850u"kg/m^3"
        A = 0.01u"m^2"
        I = 8.33e-6u"m^4"
        
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L/2, 0.0u"m", 0.0u"m"], :free)
        n3 = Node([L, 0.0u"m", 0.0u"m"], :free)
        
        sec = Section(A, E, 77u"GPa", I, I/2, 1e-6u"m^4", ρ)
        
        e1 = Element(n1, n2, sec)
        e2 = Element(n2, n3, sec)
        e1.Ψ = 0.0
        e2.Ψ = 0.0
        
        @testset "modal! auto-processes" begin
            model = FrameModel([n1, n2, n3], [e1, e2], AbstractLoad[])
            planarize!(model)
            
            # Model not processed yet
            @test model.processed == false
            
            # modal! should auto-process
            result = modal!(model; n=3)
            
            @test model.processed == true
            @test length(result.frequencies) == 3
            @test all(result.frequencies .> 0)
        end
        
        @testset "modal requires processed model" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([L/2, 0.0u"m", 0.0u"m"], :free)
            n3 = Node([L, 0.0u"m", 0.0u"m"], :free)
            
            e1 = Element(n1, n2, sec)
            e2 = Element(n2, n3, sec)
            e1.Ψ = 0.0
            e2.Ψ = 0.0
            
            model = FrameModel([n1, n2, n3], [e1, e2], AbstractLoad[])
            planarize!(model)
            
            # modal (without !) requires processing first
            @test_throws ErrorException modal(model; n=3)
            
            # After processing, it works
            process!(model)
            result = modal(model; n=3)
            @test length(result.frequencies) == 3
        end
        
        @testset "modal! and modal_analysis give same results" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([L/2, 0.0u"m", 0.0u"m"], :free)
            n3 = Node([L, 0.0u"m", 0.0u"m"], :free)
            
            e1 = Element(n1, n2, sec)
            e2 = Element(n2, n3, sec)
            e1.Ψ = 0.0
            e2.Ψ = 0.0
            
            model = FrameModel([n1, n2, n3], [e1, e2], AbstractLoad[])
            planarize!(model)
            process!(model)
            
            result1 = modal!(model; n=4)
            result2 = modal_analysis(model; n_modes=4)
            
            @test result1.frequencies ≈ result2.frequencies
            @test result1.periods ≈ result2.periods
            @test result1.mode_shapes ≈ result2.mode_shapes
        end
        
        @testset "modal! with different mass types" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([L/2, 0.0u"m", 0.0u"m"], :free)
            n3 = Node([L, 0.0u"m", 0.0u"m"], :free)
            
            e1 = Element(n1, n2, sec)
            e2 = Element(n2, n3, sec)
            e1.Ψ = 0.0
            e2.Ψ = 0.0
            
            model = FrameModel([n1, n2, n3], [e1, e2], AbstractLoad[])
            planarize!(model)
            
            result_consistent = modal!(model; n=3, mass_type=MASS_CONSISTENT)
            result_lumped = modal!(model; n=3, mass_type=MASS_LUMPED)
            
            # Both should give positive frequencies
            @test all(result_consistent.frequencies .> 0)
            @test all(result_lumped.frequencies .> 0)
            
            # Results should be similar but not identical
            @test !isapprox(result_consistent.frequencies, result_lumped.frequencies)
        end
        
        @testset "modal! with ShellModel" begin
            # Simple shell structure
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3 = Node([1.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            n4 = Node([0.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            n5 = Node([0.5u"m", 0.5u"m", 0.0u"m"], :free)  # Interior node
            
            mat = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3")
            
            shells = [
                ShellTri3((n1, n2, n5), 0.15u"m", mat),
                ShellTri3((n2, n3, n5), 0.15u"m", mat),
                ShellTri3((n3, n4, n5), 0.15u"m", mat),
                ShellTri3((n4, n1, n5), 0.15u"m", mat),
            ]
            
            model = ShellModel([n1, n2, n3, n4, n5], shells, AbstractLoad[])
            
            result = modal!(model; n=3)
            
            @test length(result.frequencies) == 3
            @test all(result.frequencies .> 0)
        end
        
        @testset "modal! with unified Model" begin
            # Mixed frame + shell model
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3 = Node([2.0u"m", 0.0u"m", 2.0u"m"], :free)
            n4 = Node([0.0u"m", 0.0u"m", 2.0u"m"], :free)
            
            # Frame columns
            sec = Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-5u"m^4", 1e-5u"m^4", 1e-6u"m^4", 7850u"kg/m^3")
            col1 = Element(n1, n4, sec)
            col2 = Element(n2, n3, sec)
            
            # Shell floor
            mat = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3")
            shell1 = ShellTri3((n1, n2, n3), 0.15u"m", mat)
            shell2 = ShellTri3((n1, n3, n4), 0.15u"m", mat)
            
            model = Model([n1, n2, n3, n4], [col1, col2], [shell1, shell2], AbstractLoad[])
            
            result = modal!(model; n=4)
            
            @test length(result.frequencies) == 4
            @test all(result.frequencies .> 0)
        end
    end
    
    @testset "ModalResult accessors" begin
        # Create simple model
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        
        sec = Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-5u"m^4", 1e-5u"m^4", 1e-6u"m^4", 7850u"kg/m^3")
        elem = Element(n1, n2, sec)
        elem.Ψ = 0.0
        
        model = FrameModel([n1, n2], [elem], AbstractLoad[])
        planarize!(model)
        
        result = modal!(model; n=3)
        
        @test length(result.frequencies) == 3
        @test length(result.periods) == 3
        @test length(result.omegas) == 3
        @test size(result.mode_shapes, 2) == 3
        @test result.n_modes == 3
        
        # Relationship between quantities
        @test result.omegas ≈ 2π .* result.frequencies
        @test result.periods ≈ 1.0 ./ result.frequencies
    end
end
