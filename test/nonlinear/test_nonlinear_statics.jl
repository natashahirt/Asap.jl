#=
Nonlinear Static Analysis Tests
===============================

Validation cases for Newton-Raphson nonlinear solver.
=#

using Test
using Asap
using Unitful

@testset "Nonlinear Static Analysis" begin
    
    @testset "Linear Case Recovery" begin
        # For a structure with small displacements, nonlinear should
        # match linear analysis
        
        L = 3.0u"m"
        E = 200e9u"Pa"
        G = 77e9u"Pa"
        
        section = Section(0.01u"m^2", E, G, 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([L, 0.0u"m", 0.0u"m"], :free)
        ]
        
        elements = [Element(nodes[1], nodes[2], section)]
        
        # Small load - linear regime
        loads = [NodeForce(nodes[2], [0.0u"N", -1e3u"N", 0.0u"N"])]
        
        # Linear analysis
        model_linear = Model(nodes, elements, loads)
        solve!(model_linear)
        u_linear = copy(model_linear.u)
        
        # Nonlinear analysis
        model_nonlin = Model(
            [Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
             Node([L, 0.0u"m", 0.0u"m"], :free)],
            [Element(nodes[1], nodes[2], section)],
            loads
        )
        result = solve_nonlinear!(model_nonlin; n_steps=5, tol=1e-6)
        
        @test result.converged
        @test length(result.load_factors) == 6  # Including initial state
        @test result.load_factors[end] ≈ 1.0
        
        # Final displacement should match linear (within tolerance)
        u_nonlin = result.displacements[end]
        @test isapprox(maximum(abs.(u_nonlin)), maximum(abs.(u_linear)), rtol=0.01)
    end
    
    @testset "Load Stepping" begin
        # Verify load stepping produces correct intermediate states
        
        L = 2.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([L, 0.0u"m", 0.0u"m"], :free)
        ]
        elements = [Element(nodes[1], nodes[2], section)]
        loads = [NodeForce(nodes[2], [0.0u"N", -10e3u"N", 0.0u"N"])]
        
        model = Model(nodes, elements, loads)
        result = solve_nonlinear!(model; n_steps=4, tol=1e-5)
        
        # Check load factors are correct increments
        @test length(result.load_factors) == 5  # 0, 0.25, 0.5, 0.75, 1.0
        @test result.load_factors[1] == 0.0
        @test isapprox(result.load_factors[end], 1.0)
        
        # Displacements should increase monotonically (linear case)
        max_disps = [maximum(abs.(u)) for u in result.displacements]
        for i in 2:length(max_disps)
            @test max_disps[i] >= max_disps[i-1]
        end
    end
    
    @testset "Convergence Tracking" begin
        L = 3.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([L, 0.0u"m", 0.0u"m"], :free)
        ]
        elements = [Element(nodes[1], nodes[2], section)]
        loads = [NodeForce(nodes[2], [0.0u"N", -5e3u"N", 0.0u"N"])]
        
        model = Model(nodes, elements, loads)
        result = solve_nonlinear!(model; n_steps=3, max_iter=10, tol=1e-6)
        
        @test result.converged
        # First step (λ=0) has 0 iterations, rest should have at least 1
        @test all(result.iterations_per_step[2:end] .> 0)
        @test all(result.equilibrium_error[2:end] .< 1e-5)  # Should converge well
    end
    
    @testset "Capacity Curve Extraction" begin
        L = 3.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([L, 0.0u"m", 0.0u"m"], :free)
        ]
        elements = [Element(nodes[1], nodes[2], section)]
        loads = [NodeForce(nodes[2], [0.0u"N", -10e3u"N", 0.0u"N"])]
        
        model = Model(nodes, elements, loads)
        result = solve_nonlinear!(model; n_steps=5)
        
        # Extract capacity curve
        disps, factors = capacity_curve(result)
        
        @test length(disps) == length(factors)
        @test factors[1] == 0.0
        @test factors[end] ≈ 1.0
        @test disps[1] == 0.0
        @test disps[end] > 0.0
    end
    
    @testset "Pushover Alias" begin
        # Test that pushover! is an alias for solve_nonlinear!
        L = 2.0u"m"
        section = Section(0.01u"m^2", 200e9u"Pa", 77e9u"Pa", 
                         1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
        
        nodes = [
            Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed),
            Node([L, 0.0u"m", 0.0u"m"], :free)
        ]
        elements = [Element(nodes[1], nodes[2], section)]
        loads = [NodeForce(nodes[2], [0.0u"N", -5e3u"N", 0.0u"N"])]
        
        model = Model(nodes, elements, loads)
        result = pushover!(model; n_steps=3)
        
        @test result isa NonlinearResult
        @test result.converged
    end
end
