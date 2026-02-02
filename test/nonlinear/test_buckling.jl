#=
Buckling Analysis Tests
=======================

Validation cases from FinEtoolsFlexStructures test_beam_buckling.jl
and classical mechanics solutions.

Reference values verified against:
- Euler buckling formula: P_cr = π²EI / (KL)²
- FinEtoolsFlexStructures benchmarks
=#

using Test
using Asap
using Unitful

@testset "Buckling Analysis" begin
    
    @testset "Cantilever Column - Euler Buckling" begin
        # Cantilever column (fixed-free): K = 2.0
        # P_cr = π²EI / (2L)² = π²EI / 4L²
        
        L = 5.0u"m"
        E = 200e9u"Pa"
        G = 77e9u"Pa"
        
        # Square section 200mm × 200mm
        b = 0.2  # m
        A = b^2 * u"m^2"
        I = b^4 / 12 * u"m^4"
        J = 0.141 * b^4 * u"m^4"  # Approx torsional constant for square
        
        section = Section(A, E, G, I, I, J)
        
        # Create cantilever column with 10 elements
        n_elem = 10
        nodes = [Node([0.0u"m", 0.0u"m", i * L / n_elem], 
                      i == 0 ? :fixed : :free) 
                 for i in 0:n_elem]
        
        elements = [Element(nodes[i], nodes[i+1], section, :column) 
                    for i in 1:n_elem]
        
        # Apply unit axial load at top (compression in -Z)
        P_applied = 1e6u"N"  # 1 MN reference load
        loads = [NodeForce(nodes[end], [0.0u"N", 0.0u"N", -P_applied])]
        
        model = Model(nodes, elements, loads)
        result = solve_buckling!(model; n=4)
        
        # Theoretical Euler buckling load
        I_val = ustrip(u"m^4", I)
        E_val = ustrip(u"Pa", E)
        L_val = ustrip(u"m", L)
        P_cr_euler = π^2 * E_val * I_val / (2.0 * L_val)^2  # K=2 for cantilever
        
        # Buckling load factor: λ × P_applied = P_cr
        λ_expected = P_cr_euler / ustrip(u"N", P_applied)
        λ_computed = result.load_factors[1]
        
        # Check within 5% (mesh refinement affects accuracy)
        @test result.n_modes >= 1
        @test isapprox(λ_computed, λ_expected, rtol=0.05)
        
        println("Cantilever buckling:")
        println("  Expected λ = $(round(λ_expected, digits=2))")
        println("  Computed λ = $(round(λ_computed, digits=2))")
    end
    
    @testset "Fixed-Pinned Column - Euler Buckling" begin
        # Fixed-pinned column: K = 0.7 (approximate)
        # P_cr = π²EI / (0.7L)²
        
        L = 4.0u"m"
        E = 200e9u"Pa"
        G = 77e9u"Pa"
        
        # Rectangular section 150mm × 100mm (weak axis buckling)
        h = 0.15  # m (height)
        b = 0.10  # m (width)
        A = h * b * u"m^2"
        Ix = b * h^3 / 12 * u"m^4"  # Strong axis
        Iy = h * b^3 / 12 * u"m^4"  # Weak axis (will buckle about this)
        J = 0.02 * (h + b) * (h * b)^3 / (h^2 + b^2) * u"m^4"  # Approx
        
        section = Section(A, E, G, Ix, Iy, J)
        
        # Create fixed-pinned column
        n_elem = 8
        nodes = Node[]
        for i in 0:n_elem
            z = i * L / n_elem
            if i == 0
                # Bottom: fixed
                push!(nodes, Node([0.0u"m", 0.0u"m", z], :fixed))
            elseif i == n_elem
                # Top: pinned (xy fixed for stability, z free for load, rotations free)
                push!(nodes, Node([0.0u"m", 0.0u"m", z], 
                             [false, false, true, true, true, true]))
            else
                push!(nodes, Node([0.0u"m", 0.0u"m", z], :free))
            end
        end
        
        elements = [Element(nodes[i], nodes[i+1], section, :column) 
                    for i in 1:n_elem]
        
        # Apply compression load
        P_applied = 500e3u"N"  # 500 kN
        loads = [NodeForce(nodes[end], [0.0u"N", 0.0u"N", -P_applied])]
        
        model = Model(nodes, elements, loads)
        result = solve_buckling!(model; n=4)
        
        # Theoretical: P_cr = π²EI_y / (KL)² with K=0.7 for fixed-pinned
        Iy_val = ustrip(u"m^4", Iy)
        E_val = ustrip(u"Pa", E)
        L_val = ustrip(u"m", L)
        K_eff = 0.7
        P_cr_euler = π^2 * E_val * Iy_val / (K_eff * L_val)^2
        
        λ_expected = P_cr_euler / ustrip(u"N", P_applied)
        λ_computed = result.load_factors[1]
        
        @test result.n_modes >= 1
        # Relax tolerance since K=0.7 is approximate
        @test isapprox(λ_computed, λ_expected, rtol=0.10)
        
        println("Fixed-pinned buckling:")
        println("  Expected λ ≈ $(round(λ_expected, digits=2)) (with K=0.7)")
        println("  Computed λ = $(round(λ_computed, digits=2))")
    end
    
    @testset "Stability Check Functions" begin
        # Simple cantilever under light load (stable)
        L = 3.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([0.0u"m", 0.0u"m", L], :free)
        ]
        
        elements = [Element(nodes[1], nodes[2], section)]
        
        # Very light load (should be very stable)
        loads = [NodeForce(nodes[2], [0.0u"N", 0.0u"N", -1000.0u"N"])]
        
        model = Model(nodes, elements, loads)
        
        @test is_stable(model)
        @test critical_load_factor(model) > 1.0
    end
    
    @testset "Result Summary" begin
        # Just test that summary functions work without error
        L = 3.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([0.0u"m", 0.0u"m", L], :free)
        ]
        elements = [Element(nodes[1], nodes[2], section)]
        loads = [NodeForce(nodes[2], [0.0u"N", 0.0u"N", -1e6u"N"])]
        
        model = Model(nodes, elements, loads)
        result = solve_buckling!(model; n=2)
        
        # Should not throw
        print_buckling_summary(result)
        @test true
    end
end
