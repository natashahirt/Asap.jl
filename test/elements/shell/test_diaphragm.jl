#=
Diaphragm Element Tests
=======================

Tests for Diaphragm elements - massless, rigid in-plane shell elements
for floor diaphragm modeling in lateral analysis.
=#

@testset "Diaphragm Elements" begin
    
    @testset "DiaphragmSection construction" begin
        # Default section
        sec = DiaphragmSection()
        @test sec.E == 1e12  # 1 TPa
        @test sec.thickness == 0.01  # 10mm
        @test sec.ν == 0.3
        @test sec.name == :diaphragm
        
        # Custom section
        sec2 = DiaphragmSection(E=200u"GPa", thickness=0.02u"m", ν=0.25, name=:steel_diaphragm)
        @test sec2.E == 200e9
        @test sec2.thickness == 0.02
        @test sec2.ν == 0.25
        @test sec2.name == :steel_diaphragm
        
        # Preset
        @test RigidDiaphragm isa DiaphragmSection
        @test RigidDiaphragm.E == 1e12
    end
    
    @testset "Diaphragm element creation" begin
        n1 = Node([0.0u"m", 0.0u"m", 3.0u"m"], :free)
        n2 = Node([6.0u"m", 0.0u"m", 3.0u"m"], :free)
        n3 = Node([6.0u"m", 6.0u"m", 3.0u"m"], :free)
        n4 = Node([0.0u"m", 6.0u"m", 3.0u"m"], :free)
        
        # Create with default section
        diaphragms = Diaphragm((n1, n2, n3, n4))
        
        @test length(diaphragms) >= 2
        @test all(d -> d isa ShellTri3, diaphragms)
        
        # Zero mass
        @test all(d -> d.ρ == 0.0, diaphragms)
        
        # High stiffness
        @test all(d -> d.E == 1e12, diaphragms)
    end
    
    @testset "Diaphragm with custom section" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Node([4.0u"m", 4.0u"m", 0.0u"m"], :free)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :free)
        
        steel_sec = DiaphragmSection(E=200u"GPa", thickness=0.005u"m")
        diaphragms = Diaphragm((n1, n2, n3, n4), steel_sec; n=3)
        
        @test all(d -> d.E == 200e9, diaphragms)
        @test all(d -> d.thickness == 0.005, diaphragms)
        @test all(d -> d.ρ == 0.0, diaphragms)  # Still massless
    end
    
    @testset "Diaphragm vs Shell mass" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        # Shell with mass
        shell_sec = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
        shells = Shell((n1, n2, n3, n4), 2, shell_sec)
        
        # Diaphragm (no mass)
        diaphragms = Diaphragm((n1, n2, n3, n4); n=2)
        
        @test all(s -> s.ρ == 2400.0, shells)
        @test all(d -> d.ρ == 0.0, diaphragms)
        
        # Total mass check
        shell_mass = sum(s.ρ * s.area * s.thickness for s in shells)
        diaphragm_mass = sum(d.ρ * d.area * d.thickness for d in diaphragms)
        
        @test shell_mass > 0
        @test diaphragm_mass == 0.0
    end
    
    @testset "Diaphragm triangular region" begin
        # Triangular diaphragm
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Node([5.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Node([2.5u"m", 4.0u"m", 0.0u"m"], :free)
        
        diaphragms = Diaphragm((n1, n2, n3); n=2)
        
        @test length(diaphragms) >= 1
        @test all(d -> d.ρ == 0.0, diaphragms)
    end
    
    @testset "Diaphragm in ShellModel" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        
        diaphragms = Diaphragm((n1, n2, n3, n4); n=2)
        all_nodes = get_nodes(diaphragms)
        
        # Find interior node for loading
        interior_nodes = filter(n -> begin
            x = ustrip(u"m", n.position[1])
            y = ustrip(u"m", n.position[2])
            x > 0.5 && x < 5.5 && y > 0.5 && y < 5.5
        end, all_nodes)
        
        if !isempty(interior_nodes)
            load_node = interior_nodes[1]
            loads = [NodeForce(load_node, [0.0u"kN", 0.0u"kN", -1.0u"kN"])]
            
            model = ShellModel(all_nodes, diaphragms, loads)
            solve!(model)
            
            # Should deflect (out-of-plane is not rigid)
            max_defl = maximum(abs(ustrip(u"m", n.displacement[3])) for n in all_nodes)
            @test max_defl > 0
        end
    end
end
