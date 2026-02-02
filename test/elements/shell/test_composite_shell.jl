using Test
using Asap
using LinearAlgebra
using Unitful

@testset "Composite Shell Elements" begin

    @testset "Ply Definition" begin
        # Test orthotropic ply
        ply = Asap.Ply("T300/5208", 181e9, 10.3e9, 7.17e9, 0.28, 0.125e-3, 0.0)
        
        @test ply.E1 == 181e9
        @test ply.E2 == 10.3e9
        @test ply.thickness == 0.125e-3
        @test ply.angle == 0.0
        
        # Check Q matrix is symmetric and positive definite
        @test ply.Qbar ≈ ply.Qbar'
        eigvals_Q = eigvals(ply.Qbar)
        @test all(eigvals_Q .> 0)
        
        # Test isotropic ply
        steel_ply = Asap.isotropic_ply("Steel", 200e9, 0.3, 0.01)
        @test steel_ply.E1 == steel_ply.E2
        @test steel_ply.Qbar[1,1] ≈ steel_ply.Qbar[2,2]
        
        @info "Ply tests passed" ply.name steel_ply.name
    end
    
    @testset "Laminate ABD Matrices" begin
        # Symmetric [0/90]s laminate - should have B ≈ 0
        ply_0 = Asap.Ply("CF", 181e9, 10.3e9, 7.17e9, 0.28, 0.125e-3, 0.0)
        ply_90 = Asap.Ply("CF", 181e9, 10.3e9, 7.17e9, 0.28, 0.125e-3, 90.0)
        
        # Symmetric layup: [0/90/90/0]
        sym_lam = Asap.Laminate("Symmetric", [ply_0, ply_90, ply_90, ply_0])
        A_sym, B_sym, D_sym = Asap.laminate_stiffnesses(sym_lam)
        
        # B should be nearly zero for symmetric laminate
        @test maximum(abs.(B_sym)) < 1e-6 * maximum(abs.(A_sym))
        
        # A and D should be symmetric and positive definite
        @test A_sym ≈ A_sym'
        @test D_sym ≈ D_sym'
        @test all(eigvals(A_sym) .> 0)
        @test all(eigvals(D_sym) .> 0)
        
        # Asymmetric layup: [0/90] should have non-zero B
        asym_lam = Asap.Laminate("Asymmetric", [ply_0, ply_90])
        A_asym, B_asym, D_asym = Asap.laminate_stiffnesses(asym_lam)
        
        # B should be non-zero for asymmetric laminate
        @test maximum(abs.(B_asym)) > 1e-6 * maximum(abs.(A_asym))
        
        @info "Laminate ABD tests passed" max_B_sym=maximum(abs.(B_sym)) max_B_asym=maximum(abs.(B_asym))
    end
    
    @testset "Laminate Thickness" begin
        ply_1 = Asap.Ply("Ply1", 100e9, 10e9, 5e9, 0.3, 0.001, 0.0)
        ply_2 = Asap.Ply("Ply2", 100e9, 10e9, 5e9, 0.3, 0.002, 45.0)
        
        lam = Asap.Laminate("Test", [ply_1, ply_2, ply_1])
        
        @test Asap.thickness(lam) ≈ 0.001 + 0.002 + 0.001
        @test lam.total_thickness ≈ 0.004
    end
    
    @testset "Transverse Shear Stiffness" begin
        ply = Asap.Ply("Test", 100e9, 10e9, 5e9, 0.3, 0.001, 0.0; G13=5e9, G23=3e9)
        lam = Asap.Laminate("Single", [ply])
        
        H = Asap.laminate_transverse_shear_stiffness(lam)
        
        @test size(H) == (2, 2)
        @test H ≈ H'  # Symmetric
        @test all(eigvals(H) .> 0)  # Positive definite
    end
    
    @testset "CompositeShellTri3 Element" begin
        # Create nodes
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)  # Equilateral
        
        # Assign node IDs
        for (i, node) in enumerate([n1, n2, n3])
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        # Create laminate
        ply = Asap.Ply("CFRP", 140e9, 10e9, 5e9, 0.3, 0.25e-3, 0.0; ρ=1600.0)
        lam = Asap.Laminate("4-ply", [ply, ply, ply, ply])  # 1mm total
        
        # Create element
        elem = Asap.CompositeShellTri3((n1, n2, n3), lam)
        Asap.process!(elem)
        
        # Check dimensions
        @test size(elem.K) == (18, 18)
        @test elem.K ≈ elem.K' rtol=1e-10  # Symmetric
        
        # Check positive semi-definiteness
        eigvals_K = eigvals(Symmetric(elem.K))
        @test all(eigvals_K .>= -1e-6 * maximum(eigvals_K))
        
        # Should have at least 6 rigid body modes (plus extras from drilling stabilization)
        rbm_count = count(x -> abs(x) < 1e-3 * maximum(abs, eigvals_K), eigvals_K)
        @test rbm_count >= 6
        
        # Check area
        # Equilateral triangle with side 1: A = √3/4
        @test elem.area ≈ sqrt(3)/4 rtol=0.01
        
        @info "CompositeShellTri3 tests passed" area=elem.area rbm_count=rbm_count
    end
    
    @testset "Isotropic vs Composite Shell Comparison" begin
        # An isotropic laminate should give same results as ShellTri3
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        E = 200e9
        ν = 0.3
        t = 0.01  # 10mm
        
        for (i, node) in enumerate([n1, n2, n3])
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        # Isotropic shell
        shell_iso = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        Asap.process!(shell_iso)
        
        # Composite shell with single isotropic ply
        ply_iso = Asap.isotropic_ply("Steel", E, ν, t)
        lam_iso = Asap.Laminate("Isotropic", [ply_iso])
        shell_comp = Asap.CompositeShellTri3((n1, n2, n3), lam_iso)
        Asap.process!(shell_comp)
        
        # Stiffness matrices should be very similar
        # (Not exactly equal due to different transverse shear formulations)
        K_iso = shell_iso.K
        K_comp = shell_comp.K
        
        # Compare membrane stiffness (DOFs 1,2 for each node)
        mem_dofs = [1, 2, 7, 8, 13, 14]
        K_iso_mem = K_iso[mem_dofs, mem_dofs]
        K_comp_mem = K_comp[mem_dofs, mem_dofs]
        
        @test K_iso_mem ≈ K_comp_mem rtol=0.1  # Within 10%
        
        @info "Isotropic comparison" max_diff=maximum(abs.(K_iso - K_comp))
    end
    
    @testset "Isotropic Ply" begin
        # Test isotropic ply creation
        steel = Asap.isotropic_ply("Steel", 200e9, 0.3, 0.01; ρ=7850.0)
        @test steel.E1 == steel.E2  # Isotropic: E1 = E2
        @test steel.thickness == 0.01
        @test steel.ρ == 7850.0
        @test steel.angle == 0.0
        
        # Q matrix should be symmetric
        @test steel.Qbar ≈ steel.Qbar'
        # And have equal diagonal terms for membrane
        @test steel.Qbar[1,1] ≈ steel.Qbar[2,2]
        
        @info "Isotropic ply test passed" name=steel.name
    end
    
    @testset "Ply Stress Recovery" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free)
        
        for (i, node) in enumerate([n1, n2, n3])
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        # Simple laminate
        ply_0 = Asap.Ply("CF", 140e9, 10e9, 5e9, 0.3, 0.5e-3, 0.0)
        ply_90 = Asap.Ply("CF", 140e9, 10e9, 5e9, 0.3, 0.5e-3, 90.0)
        lam = Asap.Laminate("Test", [ply_0, ply_90])
        
        elem = Asap.CompositeShellTri3((n1, n2, n3), lam)
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        
        # Apply membrane strain (x-direction extension)
        u = zeros(18)
        u[1] = 0.0     # node 1, u
        u[7] = 0.001   # node 2, u (1mm displacement)
        u[13] = 0.0    # node 3, u
        
        # Get ply stresses
        σ_ply1 = Asap.ply_stresses(elem, u, 1)  # Bottom ply (0°)
        σ_ply2 = Asap.ply_stresses(elem, u, 2)  # Top ply (90°)
        
        # In the 0° ply, fiber direction is aligned with load → higher σ1
        # In the 90° ply, fiber direction is perpendicular → higher σ2
        @test σ_ply1[1] > σ_ply1[2]  # 0° ply: σ1 > σ2
        @test σ_ply2[2] > σ_ply2[1]  # 90° ply: σ2 > σ1 (fibers transverse)
        
        @info "Ply stress recovery" σ_ply1=σ_ply1 σ_ply2=σ_ply2
    end
    
    @testset "AreaLoad on Composite Shell" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        for (i, node) in enumerate([n1, n2, n3])
            node.nodeID = i
            node.globalID = collect((i-1)*6+1 : i*6)
        end
        
        # Define ply directly (no presets in ASAP)
        ply = Asap.Ply("CFRP", 140e9, 10e9, 5e9, 0.3, 0.25e-3, 0.0; ρ=1600.0)
        lam = Asap.Laminate("CFRP_4ply", [ply, ply, ply, ply])
        
        elem = Asap.CompositeShellTri3((n1, n2, n3), lam)
        Asap.process!(elem)
        
        # Apply surface load
        load = Asap.AreaLoad(elem, 5000.0u"Pa"; direction=(0.0, 0.0, -1.0))
        
        # Check nodal forces
        nf = Asap.nodal_forces(load)
        total_force = sum(f[2][3] for f in nf)
        expected_force = -5000.0 * elem.area
        
        @test total_force ≈ expected_force rtol=1e-10
        @test length(nf) == 3
        
        @info "Surface load on composite shell OK" area=elem.area total_force=total_force
    end

end
