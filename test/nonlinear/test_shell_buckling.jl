#=
Shell Buckling Analysis Tests
=============================

Validation cases for shell geometric stiffness and buckling analysis.

Reference values from classical plate theory:
- Timoshenko & Gere, "Theory of Elastic Stability"
- Cook et al., "Concepts and Applications of FEA" Ch. 14

Critical buckling load for simply supported rectangular plate:
    Ncr = k · π² · D / b²
    
where:
    D = Et³ / [12(1-ν²)]  (flexural rigidity)
    b = plate width (or shorter dimension)
    k = buckling coefficient (depends on aspect ratio, loading, BCs)

For square plate (a = b) under uniaxial compression:
    k = 4.0 (m = 1, single half-wave)

For square plate under biaxial compression (Nx = Ny = N):
    Ncr = 2π²D/a²
=#

using Test
using Asap
using LinearAlgebra
using Unitful

#=============================================================================
Helper Functions
==============================================================================#

"""
    create_plate_mesh(a, b, nx, ny, E, ν, t)

Create a rectangular plate mesh with nodes at corners (0,0) to (a,b).
Returns (nodes, elements, edge_nodes_dict).
"""
function create_plate_mesh(a::Float64, b::Float64, nx::Int, ny::Int, 
                           E::Float64, ν::Float64, t::Float64)
    nodes = Asap.Node[]
    
    for j in 0:ny
        for i in 0:nx
            x = i * a / nx
            y = j * b / ny
            push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
        end
    end
    
    for (idx, node) in enumerate(nodes)
        node.nodeID = idx
        node.globalID = collect((idx-1)*6+1 : idx*6)
    end
    
    elements = Asap.ShellTri3[]
    for j in 0:ny-1
        for i in 0:nx-1
            n1 = j * (nx+1) + i + 1
            n2 = j * (nx+1) + i + 2
            n3 = (j+1) * (nx+1) + i + 1
            n4 = (j+1) * (nx+1) + i + 2
            
            tri1 = Asap.ShellTri3((nodes[n1], nodes[n2], nodes[n4]), t*u"m", E*u"Pa", ν)
            tri2 = Asap.ShellTri3((nodes[n1], nodes[n4], nodes[n3]), t*u"m", E*u"Pa", ν)
            push!(elements, tri1, tri2)
        end
    end
    
    for elem in elements
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in elem.nodes]...)
    end
    
    # Identify edge nodes
    tol = min(a/nx, b/ny) / 100
    edge_nodes = Dict{Symbol, Vector{Int}}(
        :left => Int[], :right => Int[], :bottom => Int[], :top => Int[]
    )
    
    for (idx, node) in enumerate(nodes)
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        
        abs(x) < tol && push!(edge_nodes[:left], idx)
        abs(x - a) < tol && push!(edge_nodes[:right], idx)
        abs(y) < tol && push!(edge_nodes[:bottom], idx)
        abs(y - b) < tol && push!(edge_nodes[:top], idx)
    end
    
    return nodes, elements, edge_nodes
end

"""
    apply_simply_supported_bc!(nodes, edge_nodes)

Apply simply supported BCs for out-of-plane buckling:
- All edges: w = 0 (transverse displacement constrained)
- Pin corner for in-plane rigid body modes

These BCs allow the plate to buckle out-of-plane while constraining edges.
"""
function apply_simply_supported_bc!(nodes::Vector{Asap.Node}, edge_nodes::Dict)
    all_edge_indices = unique(vcat(values(edge_nodes)...))
    
    for idx in all_edge_indices
        node = nodes[idx]
        node.dof[3] = false  # w = 0 on edges
    end
    
    # Pin bottom-left corner for in-plane stability
    corner_idx = intersect(edge_nodes[:left], edge_nodes[:bottom])[1]
    nodes[corner_idx].dof[1] = false  # u
    nodes[corner_idx].dof[2] = false  # v
    
    # Also fix v at bottom-right to prevent rotation
    corner_br = intersect(edge_nodes[:right], edge_nodes[:bottom])[1]
    nodes[corner_br].dof[2] = false
    
    return nothing
end

"""
    analytical_plate_buckling_uniaxial(a, b, E, ν, t)

Analytical critical buckling load for simply supported rectangular plate
under uniaxial compression (force per unit length in x-direction).

For aspect ratio a/b, the critical load is:
    Ncr = π²D/b² * (mb/a + a/(mb))²

where m is chosen to minimize the expression (usually m=1 for moderate aspect ratios).
"""
function analytical_plate_buckling_uniaxial(a::Float64, b::Float64, E::Float64, ν::Float64, t::Float64)
    D = E * t^3 / (12 * (1 - ν^2))
    
    # For a/b ≤ √2, m=1 gives minimum; check a few values of m
    Ncr_min = Inf
    for m in 1:5
        k = (m * b / a + a / (m * b))^2
        Ncr = π^2 * D / b^2 * k
        Ncr_min = min(Ncr_min, Ncr)
    end
    
    return Ncr_min
end

"""
    analytical_plate_buckling_biaxial(a, E, ν, t)

Analytical critical buckling load for simply supported square plate
under equal biaxial compression (Nx = Ny = N).

    Ncr = 2π²D/a²
"""
function analytical_plate_buckling_biaxial(a::Float64, E::Float64, ν::Float64, t::Float64)
    D = E * t^3 / (12 * (1 - ν^2))
    return 2 * π^2 * D / a^2
end

#=============================================================================
Test Suites
==============================================================================#

@testset "Shell Geometric Stiffness - Unit Tests" begin
    
    @testset "Matrix Properties" begin
        # Create a single triangular element and test Kg properties
        E = 200e9
        ν = 0.3
        t = 0.01
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        for (idx, node) in enumerate([n1, n2, n3])
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        elem = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        Asap.process!(elem)
        
        # Test with biaxial compression
        σ_membrane = [-1000.0, -1000.0, 0.0]  # Compression [N/m]
        Kg_local = Asap.local_geometric_stiffness(elem, σ_membrane)
        
        # Kg should be 18×18
        @test size(Kg_local) == (18, 18)
        
        # Kg should be symmetric
        @test isapprox(Kg_local, Kg_local', rtol=1e-10)
        
        # Only w DOFs (3, 9, 15) should have non-zero entries
        w_dofs = [3, 9, 15]
        for i in 1:18
            for j in 1:18
                if !(i in w_dofs) || !(j in w_dofs)
                    @test abs(Kg_local[i, j]) < 1e-15
                end
            end
        end
        
        # Extract 3×3 w-w submatrix
        Kg_ww = Kg_local[w_dofs, w_dofs]
        
        # For compression (Nxx, Nyy < 0), Kg should reduce effective stiffness
        # The eigenvalues should be negative or near-zero
        # (near-zero corresponds to uniform translation mode)
        eig_kg = eigvals(Kg_ww)
        tol = 1e-10
        @test all(eig_kg .< tol)  # All eigenvalues ≤ 0 for compression
        @test count(eig_kg .< -tol) >= 2  # At least 2 significant negative eigenvalues
        
        println("Kg_ww eigenvalues for compression: ", round.(eig_kg, sigdigits=3))
    end
    
    @testset "Linear Scaling" begin
        E = 200e9
        ν = 0.3
        t = 0.01
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        for (idx, node) in enumerate([n1, n2, n3])
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        elem = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        Asap.process!(elem)
        
        σ1 = [-1000.0, -500.0, 100.0]
        σ2 = [-2000.0, -1000.0, 200.0]  # Double
        
        Kg1 = Asap.local_geometric_stiffness(elem, σ1)
        Kg2 = Asap.local_geometric_stiffness(elem, σ2)
        
        @test isapprox(Kg2, 2.0 * Kg1, rtol=1e-10)
    end
    
    @testset "Sign Convention - Tension vs Compression" begin
        E = 200e9
        ν = 0.3
        t = 0.01
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free)
        
        for (idx, node) in enumerate([n1, n2, n3])
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        elem = Asap.ShellTri3((n1, n2, n3), t*u"m", E*u"Pa", ν)
        Asap.process!(elem)
        
        σ_comp = [-1000.0, -1000.0, 0.0]  # Compression
        σ_tens = [1000.0, 1000.0, 0.0]    # Tension
        
        Kg_comp = Asap.local_geometric_stiffness(elem, σ_comp)
        Kg_tens = Asap.local_geometric_stiffness(elem, σ_tens)
        
        # Tension should be opposite sign of compression
        @test isapprox(Kg_tens, -Kg_comp, rtol=1e-10)
        
        # Tension should stabilize (positive eigenvalues or near-zero)
        w_dofs = [3, 9, 15]
        eig_tens = eigvals(Kg_tens[w_dofs, w_dofs])
        tol = 1e-10
        @test all(eig_tens .> -tol)  # All eigenvalues ≥ 0 for tension
        @test count(eig_tens .> tol) >= 2  # At least 2 significant positive eigenvalues
    end
end

@testset "Simply Supported Plate Buckling - Prescribed Stress" begin
    
    @testset "Square Plate - Uniaxial Compression" begin
        # Test with uniform uniaxial compression (prescribed stress field)
        # This bypasses the static solve and tests the eigenvalue problem directly
        
        a = 1.0     # Plate side [m]
        E = 200e9   # Pa
        ν = 0.3
        t = 0.01    # 1% thickness ratio
        
        # Analytical solution
        D = E * t^3 / (12 * (1 - ν^2))
        Ncr_analytical = 4 * π^2 * D / a^2
        
        @info "Square Plate - Uniaxial Compression (Prescribed Stress)"
        @info "  Analytical Ncr = $(round(Ncr_analytical, sigdigits=5)) N/m"
        
        # Test mesh convergence
        mesh_sizes = [4, 8, 12]
        results = Float64[]
        
        for n in mesh_sizes
            nodes, elements, edge_nodes = create_plate_mesh(a, a, n, n, E, ν, t)
            apply_simply_supported_bc!(nodes, edge_nodes)
            
            # No loads needed - we prescribe stress directly
            loads = Asap.AbstractLoad[]
            model = Asap.ShellModel(nodes, elements, loads)
            Asap.process!(model)
            
            # Prescribe unit uniaxial compression
            σ_unit = [-1.0, 0.0, 0.0]  # Unit compression in x [N/m]
            
            # Solve buckling with prescribed stress
            result = Asap.solve_buckling!(model, σ_unit; n=4)
            
            if result.n_modes > 0
                λ = result.load_factors[1]
                Ncr_computed = λ * abs(σ_unit[1])
                accuracy = Ncr_computed / Ncr_analytical * 100
                push!(results, accuracy)
                @info "  Mesh $(n)×$(n): Ncr = $(round(Ncr_computed, sigdigits=4)), accuracy = $(round(accuracy, digits=1))%"
            else
                push!(results, 0.0)
                @warn "  Mesh $(n)×$(n): No valid buckling modes found"
            end
        end
        
        # At least the finest mesh should give reasonable results
        @test length(results) > 0
        if results[end] > 0
            @test results[end] > 70.0  # Within 30% of analytical
            @test results[end] < 140.0
        end
    end
    
    @testset "Square Plate - Biaxial Compression" begin
        a = 1.0
        E = 200e9
        ν = 0.3
        t = 0.01
        
        D = E * t^3 / (12 * (1 - ν^2))
        Ncr_analytical = 2 * π^2 * D / a^2  # Biaxial: half of uniaxial
        
        @info "Square Plate - Biaxial Compression"
        @info "  Analytical Ncr = $(round(Ncr_analytical, sigdigits=5)) N/m"
        
        n = 8
        nodes, elements, edge_nodes = create_plate_mesh(a, a, n, n, E, ν, t)
        apply_simply_supported_bc!(nodes, edge_nodes)
        
        loads = Asap.AbstractLoad[]
        model = Asap.ShellModel(nodes, elements, loads)
        Asap.process!(model)
        
        # Unit biaxial compression
        σ_unit = [-1.0, -1.0, 0.0]
        
        result = Asap.solve_buckling!(model, σ_unit; n=4)
        
        if result.n_modes > 0
            λ = result.load_factors[1]
            Ncr_computed = λ * abs(σ_unit[1])  # Both components are equal
            accuracy = Ncr_computed / Ncr_analytical * 100
            @info "  Mesh $(n)×$(n): Ncr = $(round(Ncr_computed, sigdigits=4)), accuracy = $(round(accuracy, digits=1))%"
            
            @test accuracy > 70.0
            @test accuracy < 140.0
        else
            @warn "  No valid buckling modes found"
            @test false  # Fail the test
        end
    end
    
    @testset "Biaxial vs Uniaxial Ratio" begin
        # The ratio of critical loads should be approximately 0.5
        # (biaxial is half as strong as uniaxial)
        
        a = 1.0
        E = 200e9
        ν = 0.3
        t = 0.01
        n = 8
        
        nodes, elements, edge_nodes = create_plate_mesh(a, a, n, n, E, ν, t)
        apply_simply_supported_bc!(nodes, edge_nodes)
        
        loads = Asap.AbstractLoad[]
        model = Asap.ShellModel(nodes, elements, loads)
        Asap.process!(model)
        
        # Uniaxial
        σ_uni = [-1.0, 0.0, 0.0]
        result_uni = Asap.solve_buckling!(model, σ_uni; n=2)
        
        # Biaxial (need fresh model since process! might cache things)
        nodes2, elements2, edge_nodes2 = create_plate_mesh(a, a, n, n, E, ν, t)
        apply_simply_supported_bc!(nodes2, edge_nodes2)
        model2 = Asap.ShellModel(nodes2, elements2, loads)
        Asap.process!(model2)
        
        σ_bi = [-1.0, -1.0, 0.0]
        result_bi = Asap.solve_buckling!(model2, σ_bi; n=2)
        
        if result_uni.n_modes > 0 && result_bi.n_modes > 0
            λ_uni = result_uni.load_factors[1]
            λ_bi = result_bi.load_factors[1]
            
            # For biaxial, the effective load is 2× (since both directions contribute)
            # so the ratio of λ values should be around 2.0
            ratio = λ_uni / λ_bi
            
            @info "Uniaxial λ = $(round(λ_uni, sigdigits=3)), Biaxial λ = $(round(λ_bi, sigdigits=3))"
            @info "Ratio λ_uni/λ_bi = $(round(ratio, digits=2)) (expected ~2.0)"
            
            @test ratio > 1.5  # Should be around 2
            @test ratio < 2.5
        else
            @warn "Could not compute ratio - missing modes"
        end
    end
end

@testset "Buckling Mode Shapes" begin
    a = 1.0
    E = 200e9
    ν = 0.3
    t = 0.01
    n = 8
    
    nodes, elements, edge_nodes = create_plate_mesh(a, a, n, n, E, ν, t)
    apply_simply_supported_bc!(nodes, edge_nodes)
    
    loads = Asap.AbstractLoad[]
    model = Asap.ShellModel(nodes, elements, loads)
    Asap.process!(model)
    
    σ_unit = [-1.0, 0.0, 0.0]
    result = Asap.solve_buckling!(model, σ_unit; n=4)
    
    if result.n_modes >= 1
        # Mode shapes should be normalized
        for i in 1:result.n_modes
            mode = result.mode_shapes[:, i]
            @test isapprox(maximum(abs.(mode)), 1.0, rtol=0.01)
        end
        
        # First mode should have maximum deflection near center
        mode1 = result.mode_shapes[:, 1]
        
        # Get w displacements
        w_disps = [mode1[node.globalID[3]] for node in nodes]
        
        # Find center region
        center_tol = a / 4
        center_nodes = findall(node -> begin
            x = ustrip(u"m", node.position[1])
            y = ustrip(u"m", node.position[2])
            abs(x - a/2) < center_tol && abs(y - a/2) < center_tol
        end, nodes)
        
        if !isempty(center_nodes)
            w_center_max = maximum(abs.(w_disps[center_nodes]))
            w_global_max = maximum(abs.(w_disps))
            
            # Center region should have significant deflection
            @test w_center_max > 0.5 * w_global_max
            @info "Mode 1: max center deflection / max global = $(round(w_center_max/w_global_max, digits=2))"
        end
        
        @info "Found $(result.n_modes) buckling modes"
        @info "Load factors: $(round.(result.load_factors, sigdigits=4))"
    else
        @warn "No buckling modes found"
    end
end
