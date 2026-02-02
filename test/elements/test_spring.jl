#=
Spring Element Tests
====================

Tests for grounded spring elements in Asap.jl.
Validates against analytical solutions and FinEtools reference behavior.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "Spring Element" begin
    
    @testset "Spring Constructors" begin
        # Create a test node
        n = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        
        # Single scalar stiffness (applies to all translational DOFs)
        s1 = Spring(n, 1000.0u"N/m")
        @test s1.stiffness[1:3] ≈ [1000.0, 1000.0, 1000.0]
        @test s1.stiffness[4:6] ≈ [0.0, 0.0, 0.0]
        @test s1.id == :spring
        
        # 3-element vector (translational only)
        s2 = Spring(n, [100.0u"N/m", 200.0u"N/m", 300.0u"N/m"])
        @test s2.stiffness[1:3] ≈ [100.0, 200.0, 300.0]
        @test s2.stiffness[4:6] ≈ [0.0, 0.0, 0.0]
        
        # Keyword constructor
        s3 = Spring(n; kx=500.0u"N/m", kz=1000.0u"N/m", krz=50.0u"N*m/rad")
        @test s3.stiffness[1] ≈ 500.0
        @test s3.stiffness[2] ≈ 0.0
        @test s3.stiffness[3] ≈ 1000.0
        @test s3.stiffness[6] ≈ 50.0
        
        # Custom ID
        s4 = Spring(n, 100.0u"N/m", :my_spring)
        @test s4.id == :my_spring
        
        # Dimensionless (assumes SI)
        s5 = Spring(n, 1e6)
        @test s5.stiffness[1:3] ≈ [1e6, 1e6, 1e6]
        
        # Stiffness matrix is diagonal
        @test s1.K ≈ diagm(s1.stiffness)
    end
    
    @testset "Spring Accessors" begin
        n = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        s = Spring(n; kx=100.0u"N/m", ky=200.0u"N/m", kz=300.0u"N/m", 
                      krx=10.0u"N*m/rad", kry=20.0u"N*m/rad", krz=30.0u"N*m/rad")
        
        @test spring_stiffness(s) ≈ [100.0, 200.0, 300.0, 10.0, 20.0, 30.0]
        @test translational_stiffness(s) ≈ [100.0, 200.0, 300.0]
        @test rotational_stiffness(s) ≈ [10.0, 20.0, 30.0]
    end
    
    @testset "Cantilever on Spring Support - Analytical" begin
        #=
        Test case: Cantilever beam with far end on a spring support.
        
               Fixed                    Spring
                |                         |
                |=========================|
                |          L              ↓ P
                                          k
        
        Analytical solution for displacement at spring:
        - Beam flexibility at tip: δ_beam = P*L³/(3*E*I) + P*L/(k_spring)
        - But with spring: total deflection = P*L³/(3*E*I) + P/k_spring
        
        Actually for this setup:
        - The spring adds vertical stiffness at the tip
        - Beam stiffness at tip (vertical): k_beam = 3*E*I/L³
        - Combined stiffness: k_total = k_beam + k_spring (parallel springs at same DOF)
        
        Wait, that's wrong. The spring is in series with the beam's flexibility at that point.
        
        Correct analysis:
        - Force P applied at tip
        - Deflection from beam flexibility: δ_beam = P*L³/(3EI)
        - Additional deflection from spring: δ_spring = P/k_spring
        - Total: δ_total = P*L³/(3EI) + P/k_spring
        
        Actually, if the spring is grounded (to ground), and we apply force P at the node,
        the spring just adds stiffness to that DOF. The system becomes:
        - Beam contributes stiffness k_beam = 3EI/L³ at tip in vertical direction
        - Spring contributes k_spring at same DOF
        - Total stiffness k_total = k_beam + k_spring
        - Deflection: δ = P / k_total = P / (3EI/L³ + k_spring)
        =#
        
        # Parameters
        L = 10.0u"m"
        E = 200.0u"GPa"  # Steel
        I = 1000.0u"cm^4"  # Moment of inertia
        P = 10.0u"kN"
        k_spring = 1000.0u"kN/m"
        
        # Convert to SI for calculation
        L_m = ustrip(u"m", L)
        E_Pa = ustrip(u"Pa", E)
        I_m4 = ustrip(u"m^4", I)
        P_N = ustrip(u"N", P)
        k_Nm = ustrip(u"N/m", k_spring)
        
        # Analytical solution
        k_beam = 3 * E_Pa * I_m4 / L_m^3  # Beam stiffness at tip
        k_total = k_beam + k_Nm
        δ_analytical = P_N / k_total
        
        # Create Asap model
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L, 0.0u"m", 0.0u"m"], :free)
        
        # Section: A doesn't matter for bending, use dummy values
        A = 100.0u"cm^2"
        G = 80.0u"GPa"
        J = 1.0u"cm^4"
        sec = Section(A, E, G, I, I, J)
        
        elem = Element(n1, n2, sec)
        elem.Ψ = 0.0  # Planar
        
        load = NodeForce(n2, [0.0u"N", -P, 0.0u"N"])
        
        model = Model([n1, n2], [elem], [load])
        planarize!(model)
        process!(model)
        
        # Add spring at tip (vertical direction only)
        spring = Spring(n2; ky=k_spring)
        add_springs!(model, spring)
        
        solve!(model)
        
        # Get vertical displacement (Y direction, DOF 2)
        δ_asap = abs(ustrip(u"m", n2.displacement[2]))
        
        # Compare (tolerance for numerical differences)
        @test isapprox(δ_asap, δ_analytical, rtol=1e-6)
    end
    
    @testset "Spring Dominates Flexible Beam" begin
        #=
        Test: Very flexible cantilever with a much stiffer spring at tip.
        The spring should significantly reduce deflection compared to beam-only case.
        
        Without spring: δ_beam = P*L³/(3EI)
        With spring: δ_total = P / (k_beam + k_spring) where k_beam = 3EI/L³
        =#
        
        L = 10.0u"m"
        E = 10.0u"GPa"   # Low E for flexible beam
        I = 100.0u"cm^4"  # Small I
        P = 1.0u"kN"
        
        # Make spring much stiffer than beam
        L_m = ustrip(u"m", L)
        E_Pa = ustrip(u"Pa", E)
        I_m4 = ustrip(u"m^4", I)
        P_N = ustrip(u"N", P)
        
        k_beam = 3 * E_Pa * I_m4 / L_m^3  # Beam stiffness at tip
        k_spring_val = 10 * k_beam  # Spring 10x stiffer than beam
        k_spring = k_spring_val * u"N/m"
        
        # Create model
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L, 0.0u"m", 0.0u"m"], :free)
        
        A = 100.0u"cm^2"
        G = 4.0u"GPa"
        J = 1.0u"cm^4"
        sec = Section(A, E, G, I, I, J)
        
        elem = Element(n1, n2, sec)
        elem.Ψ = 0.0
        
        load = NodeForce(n2, [0.0u"N", -P, 0.0u"N"])
        
        # Without spring
        model_no_spring = Model([n1, n2], [elem], [load])
        planarize!(model_no_spring)
        solve!(model_no_spring)
        δ_no_spring = abs(ustrip(u"m", n2.displacement[2]))
        
        # With spring - reset displacements and create new model
        n1_ws = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2_ws = Node([L, 0.0u"m", 0.0u"m"], :free)
        elem_ws = Element(n1_ws, n2_ws, sec)
        elem_ws.Ψ = 0.0
        load_ws = NodeForce(n2_ws, [0.0u"N", -P, 0.0u"N"])
        
        model_with_spring = Model([n1_ws, n2_ws], [elem_ws], [load_ws])
        planarize!(model_with_spring)
        process!(model_with_spring)
        
        spring = Spring(n2_ws; ky=k_spring)
        add_springs!(model_with_spring, spring)
        solve!(model_with_spring)
        
        δ_with_spring = abs(ustrip(u"m", n2_ws.displacement[2]))
        
        # Expected with spring
        k_total = k_beam + k_spring_val
        δ_expected = P_N / k_total
        
        # Spring should significantly reduce deflection (by ~10x since spring is 10x stiffer)
        @test isapprox(δ_with_spring, δ_expected, rtol=1e-6)
        @test δ_with_spring < δ_no_spring / 5  # Should be at least 5x less
    end
    
    @testset "Multiple Springs on Same Node" begin
        #=
        Test multiple springs on the same node.
        Springs in parallel: k_total = k1 + k2 + k_beam
        
        Use a cantilever with two springs at tip.
        =#
        
        L = 5.0u"m"
        E = 200.0u"GPa"
        I = 1000.0u"cm^4"
        P = 10.0u"kN"
        
        # Convert
        L_m = ustrip(u"m", L)
        E_Pa = ustrip(u"Pa", E)
        I_m4 = ustrip(u"m^4", I)
        P_N = ustrip(u"N", P)
        
        k_beam = 3 * E_Pa * I_m4 / L_m^3
        
        # Two springs
        k1 = 500.0u"kN/m"
        k2 = 300.0u"kN/m"
        k1_Nm = ustrip(u"N/m", k1)
        k2_Nm = ustrip(u"N/m", k2)
        
        k_total = k_beam + k1_Nm + k2_Nm
        δ_expected = P_N / k_total
        
        # Create model
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L, 0.0u"m", 0.0u"m"], :free)
        
        A = 100.0u"cm^2"
        G = 80.0u"GPa"
        J = 1.0u"cm^4"
        sec = Section(A, E, G, I, I, J)
        
        elem = Element(n1, n2, sec)
        elem.Ψ = 0.0
        
        load = NodeForce(n2, [0.0u"N", -P, 0.0u"N"])
        
        model = Model([n1, n2], [elem], [load])
        planarize!(model)
        process!(model)
        
        # Add two springs
        spring1 = Spring(n2; ky=k1)
        spring2 = Spring(n2; ky=k2)
        add_springs!(model, [spring1, spring2])
        
        solve!(model)
        
        δ_y = abs(ustrip(u"m", n2.displacement[2]))
        
        @test isapprox(δ_y, δ_expected, rtol=1e-6)
    end
    
    @testset "Unit Conversion" begin
        n = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        
        # Test various unit inputs
        s1 = Spring(n, 1.0u"kN/m")
        @test s1.stiffness[1] ≈ 1000.0  # 1 kN/m = 1000 N/m
        
        s2 = Spring(n, 1.0u"lbf/inch")
        @test isapprox(s2.stiffness[1], 175.12683, rtol=1e-4)  # 1 lbf/in ≈ 175.13 N/m
        
        s3 = Spring(n, 1.0u"kip/ft")
        @test isapprox(s3.stiffness[1], 14593.9, rtol=1e-3)  # 1 kip/ft ≈ 14594 N/m
    end
    
    @testset "Rotational Spring" begin
        #=
        Test rotational spring: cantilever with rotational spring at tip.
        The rotational spring adds moment resistance.
        =#
        
        L = 5.0u"m"
        E = 200.0u"GPa"
        I = 500.0u"cm^4"
        M = 10.0u"kN*m"  # Applied moment at tip
        k_rot = 1000.0u"kN*m/rad"  # Rotational spring stiffness
        
        # Convert for analytical
        L_m = ustrip(u"m", L)
        E_Pa = ustrip(u"Pa", E)
        I_m4 = ustrip(u"m^4", I)
        M_Nm = ustrip(u"N*m", M)
        k_rot_Nm = ustrip(u"N*m/rad", k_rot)
        
        # Beam rotational flexibility at tip under moment: θ = M*L/(E*I)
        # Beam rotational stiffness at tip: k_beam_rot = E*I/L
        k_beam_rot = E_Pa * I_m4 / L_m
        k_total_rot = k_beam_rot + k_rot_Nm
        θ_analytical = M_Nm / k_total_rot
        
        # Create model
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Node([L, 0.0u"m", 0.0u"m"], :free)
        
        A = 100.0u"cm^2"
        G = 80.0u"GPa"
        J = 1.0u"cm^4"
        sec = Section(A, E, G, I, I, J)
        
        elem = Element(n1, n2, sec)
        elem.Ψ = 0.0
        
        # Apply moment about Z axis
        load = NodeMoment(n2, [0.0u"N*m", 0.0u"N*m", M])
        
        model = Model([n1, n2], [elem], [load])
        planarize!(model)
        process!(model)
        
        # Add rotational spring about Z
        spring = Spring(n2; krz=k_rot)
        add_springs!(model, spring)
        
        solve!(model)
        
        # Get rotation about Z (DOF 6)
        θ_asap = abs(ustrip(n2.displacement[6]))  # radians
        
        @test isapprox(θ_asap, θ_analytical, rtol=1e-6)
    end
    
end
