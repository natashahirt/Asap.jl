#=
Shell Creation Tests
====================

Tests for ShellSection, Shell(), mesh(), and get_nodes().
=#

using Test
using Asap
using Asap: DT  # Access DelaunayTriangulation through Asap
using Unitful
using LinearAlgebra

@testset "Shell Creation" begin
    
    @testset "ShellSection" begin
        @testset "Direct construction" begin
            sec = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete)
            @test sec.thickness ≈ 0.15
            @test sec.E ≈ 30e9
            @test sec.ν ≈ 0.2
            @test sec.ρ ≈ 2400.0
            @test sec.name == :concrete
        end
        
        @testset "From ShellMaterial" begin
            mat = Asap.ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
            sec = ShellSection(0.2u"m", mat)
            @test sec.thickness ≈ 0.2
            @test sec.E ≈ 30e9
        end
        
        @testset "Inline Section Creation" begin
            # Presets removed - create sections inline
            sec_150 = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete_150mm)
            sec_200 = ShellSection(0.20u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete_200mm)
            @test sec_150.thickness ≈ 0.15
            @test sec_200.thickness ≈ 0.20
        end
    end
    
    @testset "mesh() - triangulation only" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Node([2.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        n4 = Node([0.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        
        # Default n=4
        tri = mesh((n1, n2, n3, n4))
        @test tri isa DT.Triangulation
        
        # Custom n
        tri2 = mesh((n1, n2, n3, n4), 6)
        @test tri2 isa DT.Triangulation
    end
    
    @testset "Shell from corners" begin
        @testset "Default n=4" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3 = Node([1.0u"m", 2.0u"m", 0.0u"m"], :fixed)
            
            sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
            shells = Shell((n1, n2, n3), sec)  # Default n=4
            
            @test length(shells) >= 1
            @test all(s.E ≈ 30e9 for s in shells)
        end
        
        @testset "Quad (auto-triangulated)" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3 = Node([4.0u"m", 3.0u"m", 0.0u"m"], :fixed)
            n4 = Node([0.0u"m", 3.0u"m", 0.0u"m"], :fixed)
            
            sec = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
            shells = Shell((n1, n2, n3, n4), sec)
            
            @test length(shells) >= 2
            @test all(s.ρ ≈ 2400.0 for s in shells)
        end
        
        @testset "Custom refinement" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3 = Node([1.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            n4 = Node([0.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            
            sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
            
            shells_coarse = Shell((n1, n2, n3, n4), 2, sec)
            
            n1b = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2b = Node([1.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3b = Node([1.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            n4b = Node([0.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            
            shells_fine = Shell((n1b, n2b, n3b, n4b), 6, sec)
            
            @test length(shells_fine) > length(shells_coarse)
        end
        
        @testset "Custom element id" begin
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :fixed)
            n3 = Node([1.0u"m", 1.0u"m", 0.0u"m"], :fixed)
            
            sec = ShellSection(0.1u"m", 30u"GPa", 0.2)
            shells = Shell((n1, n2, n3), sec; id=:floor_slab)
            
            @test all(s.id == :floor_slab for s in shells)
        end
    end
    
    @testset "Shell from external triangulation" begin
        points = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
        tri = DT.triangulate(points)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        shells = Shell(tri, sec; z=3.0)
        
        @test length(shells) >= 1
        
        # Check z-coordinate of nodes
        nodes = get_nodes(shells)
        for node in nodes
            @test ustrip(u"m", node.position[3]) ≈ 3.0
        end
    end
    
    @testset "get_nodes" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Node([2.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        n4 = Node([0.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        shells = Shell((n1, n2, n3, n4), 3, sec)
        
        nodes = get_nodes(shells)
        
        # Should include corner nodes
        @test n1 in nodes
        @test n2 in nodes
        @test n3 in nodes
        @test n4 in nodes
        
        # Should have interior nodes too (for n=3)
        @test length(nodes) > 4
        
        # No duplicates
        @test length(nodes) == length(unique(nodes))
    end
    
    @testset "Integration with ShellModel" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([2.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Node([2.0u"m", 2.0u"m", 0.0u"m"], :fixed)
        n4 = Node([0.0u"m", 2.0u"m", 0.0u"m"], :free)
        
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
        shells = Shell((n1, n2, n3, n4), sec)
        
        all_nodes = get_nodes(shells)
        loads = [AreaLoad(shells, 5.0u"kPa")]
        
        model = ShellModel(all_nodes, shells, loads)
        solve!(model)
        
        # The free corner should deflect
        @test abs(ustrip(u"m", n4.displacement[3])) > 0
    end
    
    @testset "Two-step workflow: mesh then Shell" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Node([1.0u"m", 1.0u"m", 0.0u"m"], :fixed)
        n4 = Node([0.0u"m", 1.0u"m", 0.0u"m"], :fixed)
        
        # Step 1: Create triangulation
        tri = mesh((n1, n2, n3, n4), 3)
        
        # Step 2: Create shells from triangulation
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        shells = Shell(tri, sec; z=0.0)
        
        @test length(shells) >= 2
    end
end
